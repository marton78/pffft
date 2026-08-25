#!/usr/bin/env python3
"""Wide per-library tables from long-format samples files (see docs plan 2026-08-25).

Reads pffft-bench-samples v2 files and pivots them to one CSV per
(algo, prec, xform, target-host):  rows = sessions (label/host/tree/datetime),
columns = sizes ascending, cells = median derived MFLOPS across reps.

Sparse cells are allowed (empty when a session lacks a size) -- this file is
for eyeballing, not storage; samples.py remains the source of truth.

Usage: report.py [--out-dir DIR] FILES_OR_DIRS...
"""
from __future__ import annotations

import argparse
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from samples import read_samples


def collect_sample_files(paths):
    """Expand FILES_OR_DIRS into sample .csv paths (dirs are searched recursively)."""
    out = []
    for p in map(Path, paths):
        if p.is_dir():
            out.extend(sorted(f for f in p.rglob("*.csv") if not f.name.startswith(".")))
        else:
            out.append(p)
    return sorted(dict.fromkeys(out))


def render(out_dir: Path, files) -> list[Path]:
    # group -> size -> session index -> [mflops, ...]
    # group = (algo, prec, xform, host); sessions keep input order.
    groups: dict[tuple, dict[int, dict[int, list[float]]]] = defaultdict(lambda: defaultdict(dict))
    provs: list[dict[str, str]] = []
    for i, path in enumerate(files):
        prov, rows = read_samples(path)
        provs.append(prov)
        host = prov.get("host", "")
        for s in rows:
            gkey = (s.algo, s.prec, s.xform, host)
            cell = groups[gkey][s.size].setdefault(i, [])
            cell.append(s.mflops)

    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    taken: dict[str, tuple] = {}
    for gkey in sorted(groups):
        by_size = groups[gkey]
        sessions = sorted({i for sizes in by_size.values() for i in sizes},
                          key=lambda i: provs[i].get("datetime", ""))
        sizes = sorted(by_size)
        lines = ["label,host,tree,datetime," + ",".join(str(n) for n in sizes)]
        for i in sessions:
            p = provs[i]
            row = [p.get("label", ""), p.get("host", ""), p.get("git_tree", ""),
                   p.get("datetime", "")]
            row += ["%.1f" % statistics.median(by_size[n][i]) if i in by_size[n]
                    else "" for n in sizes]
            lines.append(",".join(row))
        name = "%s-%s-%s-%s.csv" % tuple(
            re.sub(r"[^A-Za-z0-9._]+", "-", x or "unknown").strip("-")
            for x in gkey)
        if name in taken:
            stem = name[:-4]
            n = 2
            while "%s-%d.csv" % (stem, n) in taken:
                n += 1
            print("report.py: %s collides between raw keys %r and %r; "
                  "writing %s-%d.csv" % (name, taken[name], gkey, stem, n),
                  file=sys.stderr)
            name = "%s-%d.csv" % (stem, n)
        taken[name] = gkey
        dest = out_dir / name
        dest.write_text("\n".join(lines) + "\n")
        written.append(dest)
    return written


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Pivot samples files to wide MFLOPS tables.")
    ap.add_argument("--out-dir", default=".", type=Path, help="output directory (default: .)")
    ap.add_argument("paths", nargs="+", help="samples .csv files and/or directories")
    args = ap.parse_args(argv)
    files = collect_sample_files(args.paths)
    if not files:
        ap.error("no .csv files found in: %s" % ", ".join(args.paths))
    written = render(args.out_dir, files)
    for w in written:
        print(w)
    return 0


if __name__ == "__main__":
    sys.exit(main())
