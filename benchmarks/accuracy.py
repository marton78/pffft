#!/usr/bin/env python3
"""The accuracy gate: pffft at two revisions against an exact FFT, on the same inputs.

    python3 benchmarks/accuracy.py REV_A REV_B [--bits] [--target T] [--sizes LIST] [--max N]
                              [--seeds 4] [--prec LIST] [--xform LIST] [-o FILE]

A revision is anything git rev-parse understands, or '.' for the working tree.
Each is built twice as a shared library: as configured by default, and with
PFFFT_SPLIT_MAX=1024, which splits every array of more than 1 KiB into blocks
(a revision from before the split transforms builds the same library twice).
bench_fft --accuracy from the working tree loads all four and gives them the
same random inputs (benchmarks/accuracy.c). The gate fails if

  - a check within one build fails: in place against out of place, ordered
    against pffft_transform + pffft_zreorder (bit for bit), or an error of one
    input above the absolute limit (2e-6 in f32, 5e-15 in f64);
  - the split build of a revision differs from its default build in one bit;
  - B is measurably less accurate than A: at one case (precision, kind,
    length) B's error over A's exceeds 1 + --case-tol, or the geometric mean
    of that ratio over the lengths of a precision and kind exceeds
    1 + --class-tol; for the forward and for the backward transform alike.

B need not give A's bits: whether it does is reported, not judged, but with
--bits, for changes that move data and compute nothing new, B must give A's
bits in every case, and no errors are measured (many times faster). The error
is relative RMS, pooled over at least 4 * 4096 points per case, against a
double-double FFT (error about 1e-31), forward on a random signal and
backward on a random spectrum. The round trip is shown, not judged: errors of
the twiddles largely cancel in it. docs/accuracy.md has the calibration
of the tolerances.
"""
from __future__ import annotations

import argparse
import math
import os
import statistics
import sys
from pathlib import Path

import results
import targets

# the docs' lengths: powers of two to 2^22 (split under every default
# PFFFT_SPLIT_MAX), 3, 5, 9 and 15 times powers of two (radix 3, 5, 6 and the
# twiddle-free radix-8 and radix-6 last passes), and bench_fft's picked lengths
# (up to six radix-5 passes)
DEFAULT_SIZES = "all,picked"
DEFAULT_MAX = 4194304
# Tolerances, about 8 times what an arithmetic change without an effect on
# accuracy gives (3c480ca against its parent on the M2: 0.6% at most per
# case, 0.7% per class); see docs/accuracy.md.
CASE_TOL = 0.05
CLASS_TOL = 0.01


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("rev_a", metavar="REV_A", help="the baseline: a git revision, or . for the working tree")
    ap.add_argument("rev_b", metavar="REV_B", help="the revision to check against it")
    ap.add_argument("--target", default="local",
                    help="local, rosetta[:arch] (x86-64 under Rosetta on an Apple Mac, -march=arch, default "
                         "haswell), android[:serial] or ios[:udid] (default: local)")
    ap.add_argument("--sizes", default=DEFAULT_SIZES, help=f"bench_fft --sizes (default: {DEFAULT_SIZES})")
    ap.add_argument("--max", type=int, default=DEFAULT_MAX, help=f"bench_fft --max (default: {DEFAULT_MAX})")
    ap.add_argument("--seeds", type=int, default=4,
                    help="inputs for at least SEEDS * 4096 points per case, at least one (default: 4)")
    ap.add_argument("--prec", help="f32, f64 or f32,f64 (default: all built)")
    ap.add_argument("--xform", help="real, cplx or real,cplx (default: both)")
    ap.add_argument("--case-tol", type=float, default=CASE_TOL,
                    help=f"largest error ratio B/A of one case (default: {CASE_TOL})")
    ap.add_argument("--class-tol", type=float, default=CLASS_TOL,
                    help=f"largest mean error ratio B/A of a precision and kind (default: {CLASS_TOL})")
    ap.add_argument("--no-split", action="store_true", help="do not build and check the split libraries")
    ap.add_argument("--bits", action="store_true",
                    help="B must give A's bits in every case; no errors against the exact FFT (fast)")
    ap.add_argument("-v", "--verbose", action="store_true", help="print every case, not only the flagged ones")
    ap.add_argument("-o", metavar="FILE", type=Path, help="save the CSV of bench_fft --accuracy to FILE")
    args = ap.parse_args()
    if args.o and not os.access(args.o.parent, os.W_OK):
        ap.error(f"cannot write to {args.o.parent}")

    target = targets.make_target(args.target)
    revs = {"A": targets.export_rev(args.rev_a), "B": targets.export_rev(args.rev_b)}
    harness = target.build_harness(targets.REPO)
    libs = {}
    for label, (commit, src) in revs.items():
        libs[label] = target.build_lib(src, label)
        if not args.no_split:
            libs[f"{label}-split"] = target.build_lib(src, f"{label}-split", ["-DPFFFT_SPLIT_MAX=1024"], "split")
    bench_args = ["--accuracy", "--seeds", str(args.seeds), "--sizes", args.sizes, "--min", "16",
                  "--max", str(args.max)]
    if args.bits:
        bench_args.append("--bits-only")
    for option in ("prec", "xform"):
        if getattr(args, option):
            bench_args += [f"--{option}", getattr(args, option)]
    text = target.run(harness, bench_args, libs)
    if args.o:
        args.o.write_text(text)
    meta, rows = read(text)

    for label, rev in (("A", args.rev_a), ("B", args.rev_b)):
        commit = revs[label][0]
        print(f"{label} = {rev} ({'working tree' if commit == '.' else commit})")
    print(f"{meta.get('cpu')}, {meta.get('os')}, {meta.get('compiler')}, pffft {meta.get('pffft_simd')}; "
          f"at least {meta.get('seeds')} x 4096 points per case\n")
    failures = report(rows, args)
    print("\nPASSED" if not failures else f"\nFAILED: {failures}")
    sys.exit(1 if failures else 0)


def read(text: str):
    """(meta, rows) of the CSV of bench_fft --accuracy."""
    meta, rows, columns = {}, [], None
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("#"):
            key, sep, value = line[1:].partition(":")
            if sep:
                meta[key.strip()] = value.strip()
        elif line and columns is None:
            columns = line.split(",")
        elif line:
            r = dict(zip(columns, line.split(",")))
            r["n"], r["inputs"] = int(r["n"]), int(r["inputs"])
            for m in ("fwd", "bwd", "rt", "fwd_max", "bwd_max", "rt_max"):
                r[m] = float(r[m]) if r[m] else None  # fwd, bwd empty with --bits-only
            rows.append(r)
    if columns is None or columns[:2] != ["engine", "prec"] or "bwd_max" not in columns:
        raise targets.Error("bench_fft --accuracy wrote no CSV of this format")
    return meta, rows


HASHES = ("hash_fwd", "hash_bwd", "hash_un")
MEASURES = ("fwd", "bwd", "rt")
GATED = ("fwd", "bwd")


def report(rows, args) -> str:
    """Prints the tables; returns what failed, or '' if all passed."""
    by = {(r["engine"], r["prec"], r["xform"], r["n"]): r for r in rows}
    failed = []

    bad = [r for r in rows if r["checks"]]
    for r in bad:
        print(f"check failed: {r['engine']} {r['prec']} {r['xform']} n={r['n']}: {r['checks']}")
    if bad:
        failed.append(f"{len(bad)} cases fail checks within a build")

    for label in ("A", "B"):
        split = f"{label}-split"
        pairs = [(r, by[(split, *k[1:])]) for k, r in by.items() if k[0] == label and (split, *k[1:]) in by]
        differ = [r for r, s in pairs if any(r[h] != s[h] for h in HASHES)]
        if pairs:
            print(f"{split} = {label} bit for bit: {len(pairs) - len(differ)} of {len(pairs)} cases")
        for r in differ[:20]:
            print(f"  differs: {r['prec']} {r['xform']} n={r['n']}")
        if differ:
            failed.append(f"{split} differs from {label} in {len(differ)} cases")

    same = [k for k in by if k[0] == "A" and ("B", *k[1:]) in by]
    identical = sum(all(by[k][h] == by[("B", *k[1:])][h] for h in HASHES) for k in same)
    print(f"B = A bit for bit: {identical} of {len(same)} cases\n")
    if getattr(args, "bits", False):
        for k in same:
            if any(by[k][h] != by[("B", *k[1:])][h] for h in HASHES):
                print(f"  differs: {k[1]} {k[2]} n={k[3]}")
        if identical < len(same):
            failed.append(f"B differs from A in {len(same) - identical} cases (--bits)")
        return "; ".join(failed)

    for prec, xform in dict.fromkeys((r["prec"], r["xform"]) for r in rows):
        keys = sorted((k for k in same if k[1:3] == (prec, xform)), key=lambda k: k[3])
        if not keys:
            continue
        body, ratios = [], {m: [] for m in MEASURES}
        for k in keys:
            a, b = by[k], by[("B", *k[1:])]
            n = k[3]
            cells, flag = [f"{n}{'' if results.is_pow2(n) else '*'}"], ""
            for m in MEASURES:
                ratio = b[m] / a[m] if a[m] > 0 and b[m] > 0 else 1.0
                ratios[m].append(ratio)
                cells += [f"{a[m]:.3g}", f"{b[m]:.3g}", f"{ratio:.3f}"]
                if m in GATED and ratio > 1 + args.case_tol:
                    flag = " <-- B less accurate"
            if flag:
                failed.append(f"{prec} {xform} n={n}: B less accurate")
            if args.verbose or flag:
                body.append(cells + [flag])
        means = {m: gmean(v) for m, v in ratios.items()}
        worst = {m: max(v) for m, v in ratios.items()}
        largest = {m: max(by[("B", *k[1:])][m + "_max"] for k in keys) for m in MEASURES}
        print(f"{prec} {xform}: {len(keys)} lengths; error ratio B/A, geometric mean (largest): "
              + ", ".join(f"{m} {means[m]:.4f} ({worst[m]:.3f})" for m in MEASURES)
              + "; B's largest error of one input: " + ", ".join(f"{m} {largest[m]:.3g}" for m in MEASURES))
        for m in GATED:
            if means[m] > 1 + args.class_tol:
                failed.append(f"{prec} {xform}: B less accurate on average ({m} {means[m]:.4f})")
        if body:
            header = ["n"]
            for m in MEASURES:
                header += [f"A {m}", f"B {m}", "B/A"]
            print(results.table("  relative RMS error over all inputs (rt: not judged)", header + [""], body), end="")
    return "; ".join(failed)


def gmean(values) -> float:
    return math.exp(statistics.fmean(math.log(v) for v in values))


if __name__ == "__main__":
    try:
        main()
    except targets.Error as e:
        sys.exit(f"accuracy.py: {e}")
