#!/usr/bin/env python3
"""Plot bench_fft results: cost per transform / (n log2 n) against n.

    python3 bench/plot.py results.csv... [--ns] [-o out.png|out.svg]

One panel per (precision, transform kind), one color per engine. Several files
(two revisions, two devices) are overlaid, one line style per file. Without -o
the figure opens in a window. Cycles are plotted where every row has them,
else nanoseconds; --ns forces nanoseconds (compare vdsp that way, see
bench/README.md). Needs matplotlib 3.7 or later.
"""
from __future__ import annotations

import argparse
import itertools
import math
import sys
import textwrap
from pathlib import Path

import results

POW2_ONLY = {"vdsp", "green"}  # engines that support powers of two only
# Built-in engines keep their color in every figure; others (pffft@LABEL) take the colors left.
KNOWN = ["pffftu", "fftpack", "fftw-est", "fftw-meas", "kissfft", "pocketfft", "green", "vdsp", "mkl"]
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6",
          "#e34948", "#7f7f7f", "#8c564b", "#17becf", "#bcbd22"]
FILE_STYLES = ["-", "--", "-.", (0, (5, 1, 1, 1, 1, 1))]
FILE_MARKS = ["^", "v", "D", "s"]  # for the sizes that are not powers of two
HEADER = ("cpu", "os", "counter", "date", "trials")
CACHES = (("L1", "l1d_bytes"), ("L2", "l2_bytes"), ("L3", "l3_bytes"))


def next_pow2(n: int) -> int:
    return 1 << (n - 1).bit_length()


def size_label(n: float) -> str:
    return f"{n / 2**20:g}M" if n >= 2**20 else f"{n / 2**10:g}k" if n >= 2**10 else f"{n:g}"


def header(meta, keys) -> str:
    return ", ".join(f"{k} {meta.get(k, '?')}" if k in ("counter", "trials") else meta.get(k, "?") for k in keys)


def series(cases, engine, prec, xform, sizes):
    """Cost per transform / (n log2 n) of one engine in one panel: three groups
    (powers of two, other sizes, padded), each as lists [n, median, min, max].
    A power-of-two-only engine is charged its cost at the next power of two
    for the other sizes ("padded"): what its user would pay."""
    groups = ([], [], [])
    for n in sizes:
        m = n if results.is_pow2(n) or engine not in POW2_ONLY else next_pow2(n)
        if (engine, prec, xform, m) in cases:
            w = n * math.log2(n)  # floating-point log2: see commit 419a7e60cd1
            groups[0 if results.is_pow2(n) else 1 if m == n else 2].append(
                [n] + [v / w for v in results.summary(cases[engine, prec, xform, m])])
    return [[[p[c] for p in g] for c in range(4)] for g in groups]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("csv", nargs="+", help="bench_fft CSV files")
    ap.add_argument("-o", metavar="FILE", help="write the figure (.png, .svg, ...) instead of showing it")
    ap.add_argument("--ns", action="store_true", help="plot nanoseconds even where cycles were counted")
    args = ap.parse_args()
    files = []
    for p in args.csv:
        try:
            files.append((Path(p).name, *results.read_csv(Path(p))))
        except (OSError, ValueError) as e:
            sys.exit(f"plot.py: {p}: {e}")
    try:
        import matplotlib
        if args.o:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch
        from matplotlib.ticker import FixedLocator, FuncFormatter, LogLocator, NullFormatter
    except ImportError:
        sys.exit("plot.py needs matplotlib 3.7 or later: pip install matplotlib, or run it with\n"
                 "  uv run --with matplotlib python3 bench/plot.py ...")

    every = [r for _, _, rows in files for r in rows]
    field = "cycles" if not args.ns and all(r["cycles"] is not None for r in every) else "ns"
    engines = list(dict.fromkeys(r["engine"] for r in every))
    color = {e: COLORS[KNOWN.index(e)] for e in engines if e in KNOWN}
    free = [c for c in COLORS if c not in color.values()]
    for i, e in enumerate(e for e in engines if e not in color and e != "pffft"):
        color[e] = free[i % len(free)]
    color["pffft"] = "#1a1a19"
    width = {e: 2.6 if e == "pffft" else 1.9 if e.startswith("pffft") else 1.2 for e in engines}
    present = {(r["prec"], r["xform"]) for r in every}
    precs = [p for p in ("f32", "f64") if any(p == q for q, _ in present)]
    xforms = [x for x in ("real", "cplx") if any(x == y for _, y in present)]
    fig, axes = plt.subplots(len(precs), len(xforms), squeeze=False, sharex=True,
                             figsize=(6 * len(xforms) + 1, 4 * len(precs) + 1.6), layout="constrained")
    trials = max(len(v) for _, _, rows in files for v in results.by_case(rows, field).values())
    # cache lines only if all files have the same caches: a line does not say whose cache it is
    same = len({tuple(m.get(key) for _, key in CACHES) for _, m, _ in files}) == 1
    caches = [(lv, int(files[0][1][key])) for lv, key in CACHES if same and files[0][1].get(key, "").isdigit()]
    padded = cached = False
    for (i, prec), (j, xform) in itertools.product(enumerate(precs), enumerate(xforms)):
        ax = axes[i][j]
        if (prec, xform) not in present:
            ax.set_visible(False)
            if i:  # the panel above lost its shared x labels to this one
                axes[i - 1][j].tick_params(labelbottom=True)
            continue
        for k, (_, meta, rows) in enumerate(files):
            cases = results.by_case(rows, field)
            sizes = sorted({r["n"] for r in rows if (r["prec"], r["xform"]) == (prec, xform)})
            style, mark = FILE_STYLES[k % len(FILE_STYLES)], FILE_MARKS[k % len(FILE_MARKS)]
            for e in engines:
                # A line through the powers of two only: the other sizes zigzag too much to connect.
                pow2, other, pad = series(cases, e, prec, xform, sizes)
                c, lw = color[e], width[e]
                ax.plot(pow2[0], pow2[1], color=c, lw=lw, ls=style, marker="o", ms=2 + lw, mew=0, zorder=2 + lw)
                ax.plot(other[0], other[1], mark, color=c, ms=2 + lw, mew=0, zorder=2 + lw)
                ax.plot(pad[0], pad[1], mark, color=c, mfc="white", ms=2 + lw, mew=1, zorder=2)
                if trials > 1:
                    ax.fill_between(pow2[0], pow2[2], pow2[3], color=c, alpha=0.15, lw=0, zorder=1 + lw)
                    ax.vlines(other[0], other[2], other[3], color=c, alpha=0.4, lw=1, zorder=1 + lw)
                padded |= bool(pad[0])
        ns = [r["n"] for r in every if (r["prec"], r["xform"]) == (prec, xform)]
        sample = (4 if prec == "f32" else 8) * (2 if xform == "cplx" else 1)
        for level, size in caches:
            n = size / (4 * sample)  # working set of one transform: 4 arrays of n samples
            if min(ns) <= n <= max(ns):
                cached = True
                ax.axvline(n, color="#888", lw=0.8, ls=":", zorder=1)
                ax.annotate(f"{level} {size_label(size)}B", (n, 1), xycoords=("data", "axes fraction"),
                            xytext=(-1, -3), textcoords="offset points", ha="right", va="top",
                            rotation=90, fontsize=7, color="#555")
        ax.set_title(f"{prec} {'complex' if xform == 'cplx' else 'real'}", fontsize=10)
        ax.set_xscale("log", base=2)
        ax.set_yscale("log")  # keeps a 20x slower engine from squashing the rest
        ax.xaxis.set_major_locator(FixedLocator([2 ** e for e in range(4, 31)]))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: size_label(v)))
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.yaxis.set_major_locator(LogLocator(subs=(1, 1.5, 2, 3, 5, 7)))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.grid(True, which="major", color="#e6e6e6", lw=0.6)
        ax.tick_params(labelsize=8)
        if i == len(precs) - 1 or (precs[i + 1], xform) not in present:
            ax.set_xlabel("n (transform length)", fontsize=9)
        if j == 0:
            ax.set_ylabel(f"{field} per transform / (n·log2 n)", fontsize=9)

    handles = [Line2D([], [], color=color[e], lw=width[e], label=e) for e in engines]
    handles += [Line2D([], [], color="#444", marker="o", mew=0, label="n a power of two"),
                Line2D([], [], color="#444", marker="^", ls="", mew=0, label="other n")]
    if padded:
        handles.append(Line2D([], [], color="#444", marker="^", mfc="white", ls="",
                              label=f"{', '.join(sorted(POW2_ONLY & set(engines)))} padded:\n"
                                    "cost at the next power of two"))
    if trials > 1:
        handles.append(Patch(color="#444", alpha=0.15, label="min..max over trials"))
    if cached:
        handles.append(Line2D([], [], color="#888", lw=0.8, ls=":",
                              label="cache size reached by\n4 arrays of n samples\n(in, out, work, twiddles)"))
    common = [k for k in HEADER if len({m.get(k) for _, m, _ in files}) == 1]
    if len(files) > 1:
        differ = [k for k in HEADER if k not in common]
        handles += [Line2D([], [], color="#444", ls=FILE_STYLES[k % len(FILE_STYLES)],
                           marker=FILE_MARKS[k % len(FILE_MARKS)],
                           label=textwrap.fill(f"{name}: {header(meta, differ)}" if differ else name, 36))
                    for k, (name, meta, _) in enumerate(files)]
    fig.legend(handles=handles, loc="outside lower center", ncol=math.ceil(len(handles) / 3),
               fontsize=8, frameon=False)
    unit = f"{field} per transform: fastest batch" + (", median over trials" if trials > 1 else "")
    fig.suptitle("\n".join(filter(None, (header(files[0][1], common), unit))), fontsize=10)
    if args.o:
        fig.savefig(args.o, dpi=150)
    else:
        plt.show()


if __name__ == "__main__":
    main()
