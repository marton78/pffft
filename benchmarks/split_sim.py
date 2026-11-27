#!/usr/bin/env python3
"""Split thresholds that were not measured, from a sweep of those that were.

    python3 benchmarks/split_sim.py SWEEP.csv [SWEEP.csv ...] --candidate NAME=SPEC [...]
                                [--base LABEL] [--min N] [--cases REF]

A SWEEP.csv is a compare_many.py (or bench_fft) CSV whose engines pffft@LABEL
run one library with different split thresholds (--tune LABEL=split=...): a
split sweep, as tune.py and the CI rounds make. A transform's speed depends on
the split threshold only through its split depth, the number of levels the
setup splits it into blocks. So the time of any threshold at one case
(precision, kind, length) is that of a measured label with the same depth
there, in the same run (the median, if several have it): no new measurement and
no model. split_depth, a
program built from the library's own sources (split_depth.c), gives the depths.

A candidate is NAME=SPEC with SPEC as for --tune (split, split_real,
split_complex; bytes with k or m, or never), applied to the base label's
tuning. Prints, per candidate and per measured label, the geometric mean of
the speed against the base label (default: the label that never splits, if
there is one) per precision and kind, over the lengths from --min (default
65536), and how many cases no measured label could stand in for. --cases REF
lists the cases where a candidate's depth differs from that of label REF,
with the speed difference.

Example (is splitting from 1 MiB on, 1 MiB arrays included, better than from
more than 1 MiB on?):

    python3 benchmarks/split_sim.py xeon-split.csv --candidate incl=split=1048575 --cases s1m
"""
from __future__ import annotations

import argparse
import math
import re
import statistics
import subprocess
import sys
from pathlib import Path

import results
import targets

NEVER = (1 << 64) - 1


def size_value(text: str) -> int:
    if text == "never":
        return NEVER
    m = re.fullmatch(r"(\d+)([kKmM]?)", text)
    if not m:
        raise ValueError(f"not a size: {text!r}")
    return int(m.group(1)) << {"": 0, "k": 10, "K": 10, "m": 20, "M": 20}[m.group(2)]


def apply_spec(spec: str, real: int, cplx: int):
    """(split_max_real, split_max_complex) after the split keys of a --tune SPEC."""
    for item in filter(None, spec.split(",")):
        key, _, value = item.partition("=")
        if key in ("split", "split_real"):
            real = size_value(value)
        if key in ("split", "split_complex"):
            cplx = size_value(value)
    return real, cplx


def label_thresholds(meta: dict, label: str):
    """A label's split thresholds: its library's default tuning ("tuning LABEL" in the CSV header), changed
    by its --tune ("tune LABEL"); None for a library without the tuning (another revision)."""
    default = meta.get(f"tuning {label}", "")
    real = re.search(r"\bsplit_real=(\w+)", default)
    cplx = re.search(r"\bsplit_complex=(\w+)", default)
    if not real or not cplx:
        return None
    return apply_spec(meta.get(f"tune {label}", ""), size_value(real.group(1)), size_value(cplx.group(1)))


def build_split_depth() -> Path:
    build = targets.CACHE / "local" / "split-depth"
    targets.cmake_build(targets.REPO, build, targets.HARNESS_ARGS, "split_depth")
    return targets.find_file(build, ("split_depth", "split_depth.exe"))


def depths(exe: Path, prec: str, ns, thresholds) -> dict:
    """{(n, 'real'|'cplx'): [depth at each threshold]}"""
    out = subprocess.run([str(exe), prec, ",".join(str(t) for t in thresholds), *map(str, ns)],
                         capture_output=True, text=True, check=True).stdout
    d = {}
    for line in out.splitlines():
        f = line.split()
        d[(int(f[0]), f[1])] = [int(x) for x in f[2:]]
    return d


def case_ratios(rows, base: str) -> dict:
    """{label: {(prec, xform, n): median over trials of t_base / t_label}}"""
    cases = {}
    for r in rows:
        cases.setdefault((r["prec"], r["xform"], r["n"]), []).append(r)
    out = {}
    for key, same in cases.items():
        field = "cycles" if all(r["cycles"] is not None for r in same) else "ns"
        by = {}
        for r in same:
            by.setdefault(r["engine"].split("@", 1)[-1], {})[r["trial"]] = r[field]
        for label, times in by.items():
            d = [by[base][t] / times[t] for t in times if t in by.get(base, {})]
            if d:
                out.setdefault(label, {})[key] = statistics.median(d)
    return out


def gm(v) -> float:
    return (math.exp(sum(map(math.log, v)) / len(v)) - 1) * 100 if v else float("nan")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("csv", nargs="+", type=Path)
    ap.add_argument("--candidate", action="append", default=[], metavar="NAME=SPEC")
    ap.add_argument("--base", help="the label to compare with (default: the one that never splits)")
    ap.add_argument("--min", type=int, default=65536, help="the shortest length (default 65536)")
    ap.add_argument("--cases", metavar="REF", help="list the cases where a candidate's depth differs from REF's")
    args = ap.parse_args()
    exe = build_split_depth()
    for path in args.csv:
        meta, rows = results.read_csv(path)
        labels = sorted({r["engine"].split("@", 1)[-1] for r in rows if "@" in r["engine"]})
        th = {label: label_thresholds(meta, label) for label in labels}
        skipped = [l for l in labels if th[l] is None]
        labels = [l for l in labels if th[l] is not None]
        if skipped:
            print(f"(not a split sweep of this library, left out: {', '.join(skipped)})")
        base = args.base or next((l for l in labels if th[l] == (NEVER, NEVER)), None)
        if base is None:
            raise SystemExit("split_sim.py: no label never splits: give --base")
        ratios = case_ratios(rows, base)
        cands = {}
        for c in args.candidate:
            name, _, spec = c.partition("=")
            cands[name] = apply_spec(spec, *th[base])
        print(f"\n== {path}: {meta.get('cpu')}; against {base}, n >= {args.min}, geometric mean in %")
        print("  " + "label or candidate".ljust(24) + "f32 real  f32 cplx  f64 real  f64 cplx       all  (no stand-in)")
        rows_out = [(l, th[l], True) for l in labels if l != base] + [(n, t, False) for n, t in cands.items()]
        precs = sorted({k[0] for r in ratios.values() for k in r})
        ns = {p: sorted({k[2] for r in ratios.values() for k in r if k[0] == p and k[2] >= args.min}) for p in precs}
        measured = [l for l in labels]
        all_th = sorted({t for l in measured for t in th[l]} | {t for _, t, _ in rows_out for t in t})
        dep = {p: depths(exe, p, ns[p], all_th) for p in precs}

        def depth(p, x, n, pair):
            t = pair[0] if x == "real" else pair[1]
            return dep[p][(n, "real" if x == "real" else "cplx")][all_th.index(t)]

        def stand_ins(p, x, n, pair):
            """the measured labels with pair's depth at this case"""
            want = depth(p, x, n, pair)
            return [l for l in measured if depth(p, x, n, th[l]) == want and (p, x, n) in ratios.get(l, {})]

        def simulated(p, x, n, pair):
            """the median over the stand-ins (labels of the same depth run the same code)"""
            ls = stand_ins(p, x, n, pair)
            return statistics.median(ratios[l][(p, x, n)] for l in ls) if ls else None

        for name, pair, is_label in rows_out:
            per, miss = {}, 0
            for p in precs:
                for x in ("real", "cplx"):
                    for n in ns[p]:
                        v = ratios.get(name, {}).get((p, x, n)) if is_label else simulated(p, x, n, pair)
                        if v is None:
                            miss += 1 if not is_label else 0
                            continue
                        per.setdefault((p, x), []).append(v)
            cells = "".join(f"{gm(per.get((p, x), [])):+10.1f}" for p in ("f32", "f64") for x in ("real", "cplx"))
            mark = "" if is_label else " *"
            print(f"  {(name + mark).ljust(24)}{cells}{gm([v for vs in per.values() for v in vs]):+10.1f}"
                  + (f"  ({miss})" if not is_label else ""))
            if args.cases and not is_label:
                ref = args.cases
                for p in precs:
                    for x in ("real", "cplx"):
                        for n in ns[p]:
                            if depth(p, x, n, pair) == depth(p, x, n, th[ref]):
                                continue
                            a, b = simulated(p, x, n, th[ref]), simulated(p, x, n, pair)
                            if a and b:
                                print(f"      {p} {x} {n}: depth {depth(p, x, n, th[ref])} -> {depth(p, x, n, pair)}"
                                      f", {(b / a - 1) * 100:+.1f}% against {ref}")
        print("  * simulated: each case takes the time of the measured labels with the same split depth (median)")


if __name__ == "__main__":
    try:
        main()
    except targets.Error as e:
        sys.exit(f"split_sim.py: {e}")
