#!/usr/bin/env python3
"""Compare the speed of pffft at two revisions, measured side by side.

    python3 benchmarks/compare.py REV_A REV_B [--target T] [--sizes LIST] [--prec LIST]
                             [--xform LIST] [--launches 3] [--trials 2] [--fftw] [-o FILE]

A revision is anything git rev-parse understands, or '.' for the working tree.
Each is built as a shared library with its own CMake; bench_fft from the working
tree loads both and times them batch by batch, in random order, in one process.
That process is launched several times, because at small n one launch can
favour one library in all its trials (code placement, which trials do not
vary); the range then shows it. "B vs A" is B's speed-up over A, cost of A
over cost of B: positive means B is faster. -o saves the measurements of all
launches as one CSV, for plot.py.
"""
from __future__ import annotations

import argparse
import os
import shutil
import statistics
import sys
from pathlib import Path

import results
import targets

# powers of two 32 to 2^20, and from bench_fft's picked sizes one with each of
# the odd parts 3, 5, 15, 9 (the radix-3 ido == 3 path), 45, 75, 125, 375 and
# 3125 (five radix-5 passes)
DEFAULT_SIZES = "pow2,96,160,480,1152,1440,2400,4000,15360,48000,100000"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("rev_a", metavar="REV_A", help="a git revision, or . for the working tree")
    ap.add_argument("rev_b", metavar="REV_B", help="the same; B vs A > 0 means B is faster")
    ap.add_argument("--target", default="local", help="local, android[:serial] or ios[:udid] (default: local)")
    ap.add_argument("--sizes", default=DEFAULT_SIZES, help=f"bench_fft --sizes (default: {DEFAULT_SIZES})")
    ap.add_argument("--prec", help="f32, f64 or f32,f64 (default: all built)")
    ap.add_argument("--xform", help="real, cplx or real,cplx (default: both)")
    ap.add_argument("--launches", type=int, default=3, help="separate runs of bench_fft (default: 3)")
    ap.add_argument("--trials", type=int, default=2, help="trials per launch (default: 2)")
    ap.add_argument("--fftw", action="store_true", help="also measure FFTW (FFTW_MEASURE plans, n <= 65536), for reference")
    ap.add_argument("-o", metavar="FILE", type=Path, help="save the CSV of all launches to FILE")
    args = ap.parse_args()
    if args.launches < 1:
        ap.error("--launches must be at least 1")
    if args.o and not os.access(args.o.parent, os.W_OK):  # fail now, not after measuring
        ap.error(f"cannot write to {args.o.parent}")

    target = targets.make_target(args.target)
    (commit_a, src_a), (commit_b, src_b) = targets.export_rev(args.rev_a), targets.export_rev(args.rev_b)
    harness = target.build_harness(targets.REPO, fftw=args.fftw)
    libs = {"A": target.build_lib(src_a, "A"), "B": target.build_lib(src_b, "B")}
    if libs["B"] == libs["A"]:
        # a copy, so that the process loads the code twice, as for two revisions
        libs["B"] = libs["A"].with_name("copy-" + libs["A"].name)
        shutil.copy2(libs["A"], libs["B"])

    engines = ["pffft@A", "pffft@B"] + (["fftw-meas"] if args.fftw else [])
    bench_args = ["--sizes", args.sizes, "--trials", str(args.trials), "--engines", ",".join(engines)]
    for option in ("prec", "xform"):
        if getattr(args, option):
            bench_args += [f"--{option}", getattr(args, option)]
    rows = []
    for launch in range(args.launches):
        meta, launch_rows = results.read_csv(target.run(harness, bench_args, libs))
        for r in launch_rows:
            r["trial"] += launch * args.trials
        rows += launch_rows
    if args.o:
        meta.update(trials=args.launches * args.trials, launches=args.launches,
                    A=f"{args.rev_a} ({commit_a[:10]})", B=f"{args.rev_b} ({commit_b[:10]})")
        results.write_csv(args.o, meta, rows)

    for label, rev, commit in (("A", args.rev_a, commit_a), ("B", args.rev_b, commit_b)):
        print(f"{label} = {rev} ({'working tree' if commit == '.' else commit})")
    if commit_a == commit_b:
        print("A and B are the same code: the differences show the noise floor")
    print(f"{meta.get('cpu')}, {meta.get('os')}, counter: {meta.get('counter')}")
    print(f"pffft_transform_ordered, B vs A: median and range over {args.launches} launches "
          f"x {args.trials} trials; positive: B is faster\n")
    print(compare_tables(rows, engines), end="")


def compare_tables(rows, engines) -> str:
    """One table per (precision, transform kind): cost of A and B, and B vs A.
    Cycles where A and B have them in every trial, else nanoseconds."""
    cases = {f: results.by_case(rows, f) for f in ("cycles", "ns")}
    blocks, notes = [], set()
    for prec, xform in dict.fromkeys((r["prec"], r["xform"]) for r in rows):
        lines = []
        for n in sorted({r["n"] for r in rows if (r["prec"], r["xform"]) == (prec, xform)}):
            field = "cycles" if all(r["cycles"] is not None for r in rows if r["engine"] in engines[:2]
                                    and (r["prec"], r["xform"], r["n"]) == (prec, xform, n)) else "ns"
            a, b, *ref = (cases[field].get((e, prec, xform, n)) for e in engines)
            cells = [results.fmt(v and statistics.median(v)) for v in (a, b)]
            if a and b and len(a) == len(b):
                diff, lo, hi = results.summary([(x / y - 1) * 100 for x, y in zip(a, b)])  # speed-up
                cells += [f"{diff:+.2f}%", f"{lo:+.2f} .. {hi:+.2f}"]
            else:
                cells += ["-", "-"]
            lines.append((n, field, cells + [results.fmt(v and statistics.median(v)) for v in ref]))
        unit = "cycles" if any(field == "cycles" for _, field, _ in lines) else "ns"
        body = []
        for n, field, cells in lines:
            mark = ("" if results.is_pow2(n) else "*") + ("" if field == unit else " ns")
            notes.update(mark.split())
            body.append([f"{n}{mark}"] + cells)
        title = f"{prec} {xform}: {results.UNITS[unit]}, fastest batch"
        blocks.append(results.table(title, ["n", "A", "B", "B vs A", "range"] + engines[2:], body))
    texts = {"*": "* not a power of two", "ns": "ns: in nanoseconds, as cycles were not counted in every trial"}
    return "\n".join(blocks) + "".join(f"\n{t}" for m, t in texts.items() if m in notes) + "\n"


if __name__ == "__main__":
    try:
        main()
    except targets.Error as e:
        sys.exit(f"compare.py: {e}")
