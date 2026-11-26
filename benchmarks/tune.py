#!/usr/bin/env python3
"""Measure the thresholds of pffft's tuning on this CPU, and print the row of
src/pffft_tuning.h that they make.

    python3 benchmarks/tune.py [REV] [--only LIST] [--target T] [--cmake ARGS] [--max N]
                           [--launches 2] [--trials 2] [--quick] [-o DIR]

REV is the revision to measure (default: the working tree), built once as a
shared library. Each threshold of pffft_tuning_t is swept on its own, every
value a --tune of that one library (bench_fft --tune), so that all values run
the same code side by side in one process:

  split     split_max_real and split_max_complex: 128 KiB to 16 MiB and never,
            at the lengths from 65536 to 4194304 points;
  prefetch  prefetch_min (x86 only): 8 to 48 KiB and never, at 1024 to 262144;
  radix6    radix6_max: never (0) to 64 KiB, at the real lengths with a
            factor 3 from 1536 to 98304.

Per threshold, kind (real, complex) and precision it prints the speed of each
value against the default this CPU gets now (geometric mean over the lengths,
median over launches x trials), the best value, and then the table row. A
value within half a point of the default is not better: the row keeps the
default where nothing beats it clearly. Run it twice; a row should rest on
both runs. -o DIR keeps the CSVs.

--quick: one launch, one trial, every other length (a few minutes).
--cmake ARGS: CMake arguments for the builds, e.g. "-DTARGET_C_ARCH=armv7-a
-DTARGET_C_EXTRA=neon_vfpv4" for NEON on a 32-bit Raspberry Pi, whose gcc
does not enable it with the default (native). --max N: no lengths above N
(for a machine with little memory: the 4194304-point complex transforms of
all the values take about 1.5 GB in double precision).
"""
from __future__ import annotations

import argparse
import math
import re
import shlex
import statistics
import sys
from pathlib import Path

import results
import targets

NEVER = (1 << 64) - 1
KIB, MIB = 1 << 10, 1 << 20

# knob -> (tuning keys, values, kinds, sizes); 0 for radix6 is "no radix-6 pass"
SPLIT_VALUES = [128 * KIB, 256 * KIB, 512 * KIB, MIB, 2 * MIB, 4 * MIB, 8 * MIB, 16 * MIB, NEVER]
PREFETCH_VALUES = [8 * KIB, 12 * KIB, 16 * KIB, 24 * KIB, 32 * KIB, 48 * KIB, NEVER]
RADIX6_VALUES = [0, 8 * KIB, 16 * KIB, 32 * KIB, 64 * KIB]


def lengths(lo: int, hi: int, odd=(1, 3, 5, 9, 15)) -> list:
    """N = odd * 2^k in [lo, hi], valid for every transform (multiples of 32)."""
    out = set()
    for o in odd:
        n = o * 32
        while n <= hi:
            if n >= lo:
                out.add(n)
            n *= 2
    return sorted(out)


def spec_value(v: int) -> str:
    if v == NEVER:
        return "never"
    if v and v % MIB == 0:
        return f"{v // MIB}m"
    if v and v % KIB == 0:
        return f"{v // KIB}k"
    return str(v)


def c_value(v: int) -> str:
    """a threshold as src/pffft_tuning.h writes it: NEVER, BYTES(8 * MIB), BYTES(24 * KIB)"""
    if v == NEVER:
        return "NEVER"
    for unit, shift in (("MIB", 20), ("KIB", 10)):
        if v >= 1 << shift and v % (1 << shift) == 0:
            return f"BYTES({v >> shift} * {unit})"
    return f"BYTES({v})"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("rev", nargs="?", default=".", help="git revision, or . for the working tree (default)")
    ap.add_argument("--only", default="split,prefetch,radix6", help="thresholds to sweep (default: all)")
    ap.add_argument("--target", default="local", help="local, android[:serial] or ios[:udid] (default: local)")
    ap.add_argument("--launches", type=int, default=2)
    ap.add_argument("--trials", type=int, default=2)
    ap.add_argument("--prec", help="f32, f64 or f32,f64 (default: all built)")
    ap.add_argument("--quick", action="store_true", help="one launch, one trial, every other length")
    ap.add_argument("--cmake", default="", help="CMake arguments for the builds (e.g. -DTARGET_C_ARCH=...)")
    ap.add_argument("--max", type=int, default=0, help="no lengths above this")
    ap.add_argument("-o", metavar="DIR", type=Path, help="write the CSV of each sweep to DIR")
    args = ap.parse_args()
    if args.quick:
        args.launches = args.trials = 1
    knobs = [k.strip() for k in args.only.split(",") if k.strip()]
    for k in knobs:
        if k not in ("split", "prefetch", "radix6"):
            ap.error(f"--only: unknown threshold '{k}'")

    target = targets.make_target(args.target)
    commit, src = targets.export_rev(args.rev)
    extra = shlex.split(args.cmake)
    harness = target.build_harness(targets.REPO, extra)
    lib = target.build_lib(src, "tune", extra, "tune" if extra else "")
    if args.o:
        args.o.mkdir(parents=True, exist_ok=True)

    meta, best = {}, {}
    for knob in knobs:
        if knob == "split":
            sweep = dict(keys=("split_real", "split_complex"), values=SPLIT_VALUES,
                         sizes=lengths(65536, 4194304), xforms="real,cplx")
        elif knob == "prefetch":
            sweep = dict(keys=("prefetch",), values=PREFETCH_VALUES, sizes=lengths(1024, 262144),
                         xforms="real,cplx")
        else:
            sweep = dict(keys=("radix6",), values=RADIX6_VALUES, sizes=lengths(1536, 98304, odd=(3, 9, 15)),
                         xforms="real")
        sizes = [n for n in sweep["sizes"] if not args.max or n <= args.max]
        sizes = sizes[::2] if args.quick else sizes
        meta, rows = run_sweep(target, harness, lib, knob, sweep, sizes, args)
        if knob == "prefetch" and not meta.get("pffft_cpu", "").startswith("x86"):
            print("prefetch: not an x86 CPU, the prefetches are compiled out: skipped\n")
            continue
        if args.o:
            results.write_csv(args.o / f"tune-{knob}.csv", meta, rows)
        best[knob] = report(knob, sweep, meta, rows)

    print_row(meta, best)


def run_sweep(target, harness, lib, knob, sweep, sizes, args):
    """One library under a label per value, plus the default; returns (meta, rows)."""
    libs = {"default": lib}
    bench_args = ["--sizes", ",".join(map(str, sizes)), "--xform", sweep["xforms"], "--trials", str(args.trials)]
    if args.prec:
        bench_args += ["--prec", args.prec]
    for v in sweep["values"]:
        label = f"{knob}-{spec_value(v)}"
        libs[label] = lib
        bench_args += ["--tune", f"{label}=" + ",".join(f"{k}={spec_value(v)}" for k in sweep["keys"])]
    rows, meta = [], {}
    for launch in range(args.launches):
        print(f"{knob}: launch {launch + 1} of {args.launches}", file=sys.stderr)
        meta, launch_rows = results.read_csv(target.run(harness, bench_args, libs))
        for r in launch_rows:
            r["trial"] += launch * args.trials
        rows += launch_rows
    meta.update(trials=args.launches * args.trials, launches=args.launches)
    return meta, rows


def default_value(meta, knob, xform) -> int | None:
    """The default this CPU gets now, from the CSV header's "tuning default"."""
    text = meta.get("tuning default", "")
    key = {"split": f"split_{'real' if xform == 'real' else 'complex'}", "prefetch": "prefetch",
           "radix6": "radix6"}[knob]
    m = re.search(rf"\b{key}=(\w+)", text)
    if not m:
        return None
    return NEVER if m.group(1) == "never" else int(m.group(1))


def gains(rows, base: str, engine: str, prec: str, xform: str) -> list:
    """Per length, the median over trials of time(base) / time(engine) - 1, in %."""
    out = []
    for n in sorted({r["n"] for r in rows if (r["prec"], r["xform"]) == (prec, xform)}):
        same = [r for r in rows if (r["prec"], r["xform"], r["n"]) == (prec, xform, n)]
        field = "cycles" if all(r["cycles"] is not None for r in same) else "ns"
        by = {}
        for r in same:
            by.setdefault(r["engine"], {})[r["trial"]] = r[field]
        d = [by[base][t] / by[engine][t] for t in by.get(base, {}) if t in by.get(engine, {})]
        if d:
            out.append(statistics.median(d))
    return out


def geomean_pct(ratios: list) -> float:
    return (math.exp(sum(map(math.log, ratios)) / len(ratios)) - 1) * 100 if ratios else float("nan")


def report(knob, sweep, meta, rows) -> dict:
    """Prints the table of a sweep; returns xform -> best value (pooled over the precisions)."""
    base = "pffft@default"
    best = {}
    for xform in sorted({r["xform"] for r in rows}):
        precs = sorted({r["prec"] for r in rows if r["xform"] == xform})
        default = default_value(meta, knob, xform)
        print(f"\n{knob}, {xform}: speed against the default ({spec_value(default) if default is not None else '?'})"
              f", geometric mean over the lengths, %; >0: faster")
        print("value".rjust(10) + "".join(p.rjust(10) for p in precs) + "all".rjust(10))
        scores = {}
        for v in sweep["values"]:
            engine = f"pffft@{knob}-{spec_value(v)}"
            per = [gains(rows, base, engine, p, xform) for p in precs]
            pooled = geomean_pct([g for gs in per for g in gs])
            scores[v] = pooled
            mark = "  (default)" if v == default else ""
            print(spec_value(v).rjust(10) + "".join(f"{geomean_pct(g):+10.2f}" for g in per) + f"{pooled:+10.2f}{mark}")
        best[xform] = choose(scores, default)
    if knob != "split" and len(best) > 1:
        # one threshold for both kinds: pooled over them
        default = default_value(meta, knob, "real")
        scores = {}
        for v in sweep["values"]:
            engine = f"pffft@{knob}-{spec_value(v)}"
            scores[v] = geomean_pct([g for xform in best for p in sorted({r["prec"] for r in rows})
                                     for g in gains(rows, base, engine, p, xform)])
        print(f"\n{knob}, both kinds:")
        best["real"] = choose(scores, default)
    return best


def choose(scores: dict, default) -> int:
    """The best value; the default unless another is better by more than half a point."""
    top = max(scores, key=lambda v: scores[v])
    pick = top if default is None or scores[top] > scores.get(default, 0.0) + 0.5 else default
    print(f"best: {spec_value(top)} ({scores[top]:+.2f}%); for the table: {spec_value(pick)}")
    return pick


def print_row(meta, best) -> None:
    cpu = meta.get("pffft_cpu", "")
    print(f"\npffft sees: {cpu or 'nothing (a revision from before the tuning?)'}")
    m_x86 = re.match(r"x86 (\w+) family 0x([0-9a-f]+) model 0x([0-9a-f]+)", cpu)
    m_arm = re.match(r"Arm implementer 0x([0-9a-f]+) part 0x([0-9a-f]+)", cpu)
    if m_x86:
        vendor = {"GenuineIntel": "X86_INTEL", "AuthenticAMD": "X86_AMD", "HygonGenuine": "X86_HYGON"}.get(
            m_x86.group(1), "0")
        ident = f"0x{int(m_x86.group(2), 16):02x}{int(m_x86.group(3), 16):02x}"
        head = f'PFFFT_CPU_X86, {vendor}, {ident}, {ident}'
    elif m_arm:
        part = f"0x{int(m_arm.group(2), 16):03x}"
        head = f'PFFFT_CPU_ARM, 0x{int(m_arm.group(1), 16):02x}, {part}, {part}'
    else:
        head = "PFFFT_CPU_OTHER, 0, 0, 0xFFFFFFFF"
    dflt = {}
    for key, knob, xform in (("split_real", "split", "real"), ("split_complex", "split", "cplx"),
                             ("prefetch", "prefetch", "real"), ("radix6", "radix6", "real")):
        v = best.get(knob, {}).get(xform)
        dflt[key] = v if v is not None else default_value(meta, knob, xform)
    vals = ", ".join(c_value(dflt[k]) if dflt[k] is not None else "?" for k in
                     ("split_real", "split_complex", "prefetch", "radix6"))
    print(f"machine: {meta.get('cpu')}, {meta.get('os')}, {meta.get('compiler')}")
    print("row for src/pffft_tuning.h (above the row that matches now; ids as found, widen them\n"
          "to the family or the core type if you know it):\n")
    print(f'  {{"{meta.get("cpu", "?")}", {head},\n   {vals}}},')


if __name__ == "__main__":
    try:
        main()
    except targets.Error as e:
        sys.exit(f"tune.py: {e}")
