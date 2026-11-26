#!/usr/bin/env python3
"""Time many pffft revisions side by side in one bench_fft process.

    python3 benchmarks/compare_many.py REV [REV ...] [--tune LABEL=SPEC ...] [--sizes LIST]
                                  [--prec LIST] [--xform LIST] [--launches 3] [--trials 2]
                                  [--target T] [-o FILE]

A REV is a git revision, optionally LABEL=REV. Each distinct commit is built once
as a shared library; a revision given twice is loaded twice (as a copy), which
measures the noise floor inside the run. The engines are pffft@LABEL (LABEL
defaults to the revision as given). Prints, per case, each revision's speed-up
over the first one (median and range over launches x trials).

--tune LABEL=SPEC makes the setups of pffft@LABEL with the library's default
tuning changed by SPEC (bench_fft --help: e.g. split=8m,prefetch=never). A
label with a --tune loads the library of its commit itself, not a copy, so
that tunings of one commit compare the same code: e.g.
"base=HEAD s8m=HEAD --tune s8m=split=8m".

--cmake LABEL="ARGS" builds the library of LABEL with these CMake arguments
added, e.g. --cmake out="-DPFFFT_PREFETCH_DIST=0"; --cmake "*=ARGS" adds
them to every label (e.g. -DTARGET_C_ARCH=x86-64-v3 for a generic build).
Labels with the same commit and arguments share one build.
"""
from __future__ import annotations

import argparse
import hashlib
import shlex
import shutil
import statistics
import sys
from pathlib import Path

import results
import targets

DEFAULT_SIZES = "1024,2048,4000,4096,8192,15360,16384,32768,48000,65536,100000,131072,262144"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("revs", nargs="+", metavar="REV")
    ap.add_argument("--tune", action="append", default=[], metavar="LABEL=SPEC",
                    help="tuning of the setups of pffft@LABEL (repeatable)")
    ap.add_argument("--cmake", action="append", default=[], metavar="LABEL=ARGS",
                    help="CMake arguments for the library of LABEL, or of every label with * (repeatable)")
    ap.add_argument("--sizes", default=DEFAULT_SIZES)
    ap.add_argument("--prec")
    ap.add_argument("--xform", default="real")
    ap.add_argument("--launches", type=int, default=3)
    ap.add_argument("--trials", type=int, default=2)
    ap.add_argument("-o", metavar="FILE", type=Path)
    ap.add_argument("--target", default="local", help="local, android[:serial] or ios[:udid] (default: local)")
    args = ap.parse_args()

    tuned = {t.partition("=")[0] for t in args.tune}
    cmake = {}
    for c in args.cmake:
        label, sep, extra = c.partition("=")
        if not sep:
            ap.error(f"--cmake: expected LABEL=ARGS, got '{c}'")
        cmake.setdefault(label, []).extend(shlex.split(extra))
    target = targets.make_target(args.target)
    harness = target.build_harness(targets.REPO)
    libs, commits, seen = {}, {}, {}
    for spec in args.revs:
        label, _, rev = spec.rpartition("=")
        label = label or rev
        if label in libs:
            label += "-copy"
        commit, src = targets.export_rev(rev)
        extra = cmake.get("*", []) + cmake.get(label, [])
        # a build directory per set of arguments; the variant names it
        variant = hashlib.sha1(" ".join(extra).encode()).hexdigest()[:10] if extra else ""
        path = target.build_lib(src, label.replace("/", "_"), extra, variant)
        if path in seen.values() and label not in tuned:  # same library: load a copy, as for another revision
            copy = path.with_name(f"copy-{label}-{path.name}")
            shutil.copy2(path, copy)
            path = copy
        seen[label] = path
        libs[label] = path
        commits[label] = commit
    engines = [f"pffft@{label}" for label in libs]
    bench_args = ["--sizes", args.sizes, "--trials", str(args.trials), "--engines", ",".join(engines)]
    for t in args.tune:
        bench_args += ["--tune", t]
    for option in ("prec", "xform"):
        if getattr(args, option):
            bench_args += [f"--{option}", getattr(args, option)]
    rows, meta = [], {}
    for launch in range(args.launches):
        meta, launch_rows = results.read_csv(target.run(harness, bench_args, libs))
        for r in launch_rows:
            r["trial"] += launch * args.trials
        rows += launch_rows
    meta.update(trials=args.launches * args.trials, launches=args.launches,
                revisions=" ".join(f"{l}={c[:10]}" for l, c in commits.items()))
    if args.tune:
        meta.update(tunes=" ".join(args.tune))
    if args.cmake:
        meta.update(cmake=" ; ".join(args.cmake))
    if args.o:
        results.write_csv(args.o, meta, rows)
    print(f"{meta.get('cpu')}, {meta.get('os')}, counter: {meta.get('counter')}")
    print(report(rows, engines))


def report(rows, engines) -> str:
    """Per (prec, xform, n): each engine's speed-up over the first, in % (median [min, max])."""
    out = []
    base = engines[0]
    for prec, xform in dict.fromkeys((r["prec"], r["xform"]) for r in rows):
        out.append(f"\n{prec} {xform}: speed-up over {base}, % (median, range); >0: faster")
        out.append("n".rjust(8) + "".join(e.removeprefix("pffft@").rjust(22) for e in engines[1:]))
        for n in sorted({r["n"] for r in rows if (r["prec"], r["xform"]) == (prec, xform)}):
            same = [r for r in rows if (r["prec"], r["xform"], r["n"]) == (prec, xform, n)]
            field = "cycles" if all(r["cycles"] is not None for r in same) else "ns"
            by = {}
            for r in same:
                by.setdefault(r["engine"], {})[r["trial"]] = r[field]
            line = str(n).rjust(8)
            for e in engines[1:]:
                d = [(by[base][t] / by[e][t] - 1) * 100 for t in by.get(base, {}) if t in by.get(e, {})]
                line += (f"{statistics.median(d):+.2f} [{min(d):+.1f},{max(d):+.1f}]" if d else "-").rjust(22)
            out.append(line + ("" if field == "cycles" else "  (ns)"))
    return "\n".join(out)


if __name__ == "__main__":
    try:
        main()
    except targets.Error as e:
        sys.exit(f"compare_many.py: {e}")
