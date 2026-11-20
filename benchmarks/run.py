#!/usr/bin/env python3
"""Build bench_fft from the working tree for a target, run it, print the results.

    python3 benchmarks/run.py [--target T] [--fftw] [-o FILE] [-- bench_fft options]

Example: python3 benchmarks/run.py -o results.csv -- --sizes 256,1024,1152 --prec f32
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import results
import targets


def main() -> None:
    argv = sys.argv[1:]
    split = argv.index("--") if "--" in argv else len(argv)
    passthrough = argv[split + 1:]
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0],
                                 usage="%(prog)s [--target T] [--fftw] [-o FILE] [-- bench_fft options]",
                                 epilog="Options after -- go to bench_fft (see its --help); run.py adds --csv -.")
    ap.add_argument("--target", default="local", help="local, android[:serial] or ios[:udid] (default: local)")
    ap.add_argument("--fftw", action="store_true", help="build FFTW for the target and measure it too")
    ap.add_argument("-o", metavar="FILE", type=Path, help="save the CSV to FILE")
    args = ap.parse_args(argv[:split])
    if args.o and not os.access(args.o.parent, os.W_OK):  # fail now, not after measuring
        ap.error(f"cannot write to {args.o.parent}")

    target = targets.make_target(args.target)
    harness = target.build_harness(targets.REPO, fftw=args.fftw)
    csv = target.run(harness, passthrough)
    if args.o:
        args.o.write_text(csv)
    meta, rows = results.read_csv(csv)
    show = passthrough[passthrough.index("--show") + 1] if "--show" in passthrough[:-1] else "cycles"
    if show == "cycles" and any(r["cycles"] is None for r in rows):  # as plot.py does
        show = "ns"
    print("".join(f"{k}: {v}\n" for k, v in meta.items() if k != "format"))
    print(results.engine_tables(rows, show), end="")


if __name__ == "__main__":
    try:
        main()
    except targets.Error as e:
        sys.exit(f"run.py: {e}")
