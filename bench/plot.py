#!/usr/bin/env python3
"""Plot bench_fft results.

    python3 bench/plot.py results.csv... [--engines LIST] [--gflops] [-o out.png|out.svg]

A panel per (precision, transform kind), with n on the x axis:

- for bench_fft CSVs, the speed of each engine, higher is faster: flops per
  CPU cycle, counting 5·n·log2(n) flops per complex transform and 2.5·n·log2(n)
  per real one, as benchFFT does. A solid line joins the powers of two, a
  dashed one the other lengths. With --gflops, or where cycles were not
  counted, GFLOPS (flops per ns), which depend on the clock;
- for a CSV that compare.py wrote (-o): B vs A in percent, the median and the
  range over all launches and trials.

Several files (two devices, say) are overlaid, one marker shape per file.
Without -o the figure opens in a window. Needs matplotlib 3.7 or later.
"""
from __future__ import annotations

import argparse
import itertools
import math
import sys
from pathlib import Path

import results

POW2_ONLY = {"vdsp", "green"}  # engines that support powers of two only
# An engine keeps its color in every figure, in the order of the palette;
# pffft@LABEL engines take the colors left over.
SERIES = {"pffft": "#2a78d6", "fftw-meas": "#eb6834", "vdsp": "#1baf7a", "pocketfft": "#eda100",
          "kissfft": "#e87ba4", "fftpack": "#008300", "green": "#6250d6", "fftw-est": "#e34948",
          "mkl": "#898781"}
SPARE = ["#e34948", "#6250d6", "#008300", "#e87ba4", "#eda100", "#1baf7a"]
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
SURFACE, GRID, AXIS = "#fcfcfb", "#e1e0d9", "#c3c2b7"
FASTER, SLOWER = "#2a78d6", "#e34948"
FILE_MARKS = ["o", "s", "D", "^"]
HEADER = ("cpu", "os", "counter")
CACHES = (("L1", "l1d_bytes"), ("L2", "l2_bytes"), ("L3", "l3_bytes"))
PANEL_TITLES = {"real": "real", "cplx": "complex"}


def size_label(n: float) -> str:
    return f"{n / 2**20:g}M" if n >= 2**20 else f"{n / 2**10:g}k" if n >= 2**10 else f"{n:g}"


def bytes_label(b: int) -> str:
    return f"{b / 2**20:g} MiB" if b >= 2**20 else f"{b / 2**10:g} KiB"


def machine(meta) -> str:
    os_name = meta.get("os", "")
    return ", ".join(x for x in (meta.get("cpu"), os_name) if x) or "unknown machine"


def style_axes(ax) -> None:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=MUTED, labelcolor=INK2, labelsize=8, length=3, width=0.6)
    ax.tick_params(which="minor", length=0)
    ax.grid(True, which="major", color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def size_axis(ax, sizes) -> None:
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
    ax.set_xscale("log", base=2)
    lo, hi = min(sizes), max(sizes)
    ax.set_xlim(lo / 2**0.3, hi * 2**0.3)
    ticks = [2**e for e in range(1, 31) if lo / 2 <= 2**e <= hi * 2]
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: size_label(v)))
    ax.xaxis.set_minor_formatter(NullFormatter())


def log_ticks(lo: float, hi: float):
    """Round values for a log axis from lo to hi: the sparsest set that
    gives at least 4 ticks."""
    sets = [[round(m * 10**d, 6) for d in range(-3, 7) for m in ms]
            for ms in ((1, 2, 5), (1, 1.5, 2, 3, 5, 7), (1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8, 9))]
    ticks = []
    for candidates in sets:
        ticks = [t for t in candidates if lo <= t <= hi]
        if len(ticks) >= 4:
            break
    return ticks


def load(paths):
    files = []
    for p in paths:
        try:
            files.append((Path(p).stem, *results.read_csv(Path(p))))
        except (OSError, ValueError) as e:
            sys.exit(f"plot.py: {p}: {e}")
    return files


def panels(every):
    present = {(r["prec"], r["xform"]) for r in every}
    precs = [p for p in ("f32", "f64") if any(p == q for q, _ in present)]
    xforms = [x for x in ("real", "cplx") if any(x == y for _, y in present)]
    return present, precs, xforms


class Layout:
    """A grid of panels under a title, a subtitle and a legend, with a key
    below: sizes in inches, so that text keeps its room at any grid size."""
    PANEL_W, PANEL_H, LEFT, RIGHT, ROW = 6.0, 3.6, 1.05, 0.3, 0.23

    def __init__(self, plt, precs, xforms, legend_rows: int, key_rows: int):
        rows, cols = len(precs), len(xforms)
        self.top = 0.95 + legend_rows * self.ROW + 0.35  # title, subtitle, legend, panel title
        self.bottom = 0.55 + key_rows * self.ROW + 0.15  # x labels, key
        gap_w, gap_h = 1.05, 0.75
        self.w = self.LEFT + cols * self.PANEL_W + (cols - 1) * gap_w + self.RIGHT
        self.h = self.top + rows * self.PANEL_H + (rows - 1) * gap_h + self.bottom
        self.fig, self.axes = plt.subplots(rows, cols, squeeze=False, sharex=True, figsize=(self.w, self.h))
        self.fig.set_facecolor(SURFACE)
        self.fig.subplots_adjust(left=self.LEFT / self.w, right=1 - self.RIGHT / self.w,
                                 top=1 - self.top / self.h, bottom=self.bottom / self.h,
                                 wspace=gap_w / self.PANEL_W, hspace=gap_h / self.PANEL_H)

    def y(self, inches_from_top: float) -> float:
        return 1 - inches_from_top / self.h

    @staticmethod
    def by_rows(handles, ncol: int):
        """matplotlib fills a legend column by column: reorder to read by rows."""
        rows = math.ceil(len(handles) / ncol)
        grid = [handles[r * ncol:(r + 1) * ncol] for r in range(rows)]
        return [grid[r][c] for c in range(ncol) for r in range(rows) if c < len(grid[r])]

    def header(self, title: str, subtitle: str, handles, ncol: int) -> None:
        x = 0.25 / self.w
        self.fig.text(x, self.y(0.22), title, ha="left", va="top", fontsize=14, fontweight="bold", color=INK)
        self.fig.text(x, self.y(0.55), subtitle, ha="left", va="top", fontsize=9, color=INK2)
        if handles:
            self.fig.legend(handles=self.by_rows(handles, ncol), loc="upper left", bbox_to_anchor=(x - 0.08 / self.w, self.y(0.82)),
                            ncol=ncol, frameon=False, fontsize=9, handlelength=2.6, columnspacing=2.0,
                            labelcolor=INK, borderaxespad=0, borderpad=0)

    def key(self, handles, ncol: int) -> None:
        self.fig.legend(handles=self.by_rows(handles, ncol), loc="lower left",
                        bbox_to_anchor=(0.17 / self.w, 0.12 / self.h),
                        ncol=ncol, frameon=False, fontsize=8, handlelength=2.4, columnspacing=2.0,
                        labelcolor=INK2, borderaxespad=0, borderpad=0)


def engine_colors(engines):
    color = {e: SERIES[e] for e in engines if e in SERIES}
    spare = [c for c in SPARE if c not in color.values()]
    for e in engines:
        if e not in color:
            color[e] = spare.pop(0) if spare else MUTED
    return color


def label_for(e: str) -> str:
    return {"fftw-meas": "FFTW (measure)", "fftw-est": "FFTW (estimate)",
            "vdsp": "vDSP (Accelerate)", "pocketfft": "PocketFFT", "kissfft": "KISS FFT",
            "fftpack": "FFTPACK", "green": "Green FFT", "mkl": "Intel MKL"}.get(e, e)


def plot_engines(args, files, plt) -> None:
    """Speed against n: flops per cycle, or GFLOPS with --gflops."""
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter, NullLocator

    every = [r for _, _, rows in files for r in rows]
    field = "cycles" if not args.gflops and all(r["cycles"] is not None for r in every) else "ns"
    engines = list(dict.fromkeys(r["engine"] for r in every))
    # pffft first in the legend, drawn last (on top)
    engines.sort(key=lambda e: (not e.startswith("pffft"), list(SERIES).index(e) if e in SERIES else 99))
    color = engine_colors(engines)
    width = {e: 2.4 if e.split("@")[0] == "pffft" else 1.5 for e in engines}
    present, precs, xforms = panels(every)
    ncol = min(len(engines), 5)
    key_ncol = 5 if len(files) == 1 else 4
    key_rows = math.ceil((5 + (len(files) if len(files) > 1 else 0)) / key_ncol)  # at most this many keys
    single = len(engines) == 1 and len(files) == 1  # no legend: the title names the engine
    layout = Layout(plt, precs, xforms, 0 if single else math.ceil(len(engines) / ncol), key_rows)
    fig, axes = layout.fig, layout.axes
    trials = max(len(v) for _, _, rows in files for v in results.by_case(rows, field).values())
    same = len({tuple(m.get(key) for _, key in CACHES) for _, m, _ in files}) == 1
    caches = [(lv, int(files[0][1][key])) for lv, key in CACHES if same and files[0][1].get(key, "").isdigit()]
    cached = other_n = False

    for (i, prec), (j, xform) in itertools.product(enumerate(precs), enumerate(xforms)):
        ax = axes[i][j]
        style_axes(ax)
        if (prec, xform) not in present:
            ax.set_visible(False)
            continue
        all_sizes = sorted({r["n"] for r in every if (r["prec"], r["xform"]) == (prec, xform)})
        ys = []
        for k, (_, meta, rows) in enumerate(files):
            cases = results.by_case(rows, field)
            sizes = sorted({r["n"] for r in rows if (r["prec"], r["xform"]) == (prec, xform)})
            mark = FILE_MARKS[k % len(FILE_MARKS)]
            for e in engines[::-1]:
                pts = ([], [])  # powers of two, other n
                for n in sizes:
                    v = cases.get((e, prec, xform, n))
                    if not v:
                        continue
                    # floating-point log2, see commit 419a7e60cd1
                    flops = (5 if xform == "cplx" else 2.5) * n * math.log2(n)
                    vals = [flops / x for x in v]  # per cycle, or per ns = GFLOPS
                    pts[0 if results.is_pow2(n) else 1].append((n, *results.summary(vals)))
                c, lw = color[e], width[e]
                for group, ls, alpha in ((0, "-", 0.14), (1, (0, (3.5, 2)), 0.08)):
                    if not pts[group]:
                        continue
                    g = list(zip(*pts[group]))
                    lwg = lw if group == 0 else max(1.1, lw * 0.6)
                    ax.plot(g[0], g[1], color=c, lw=lwg, ls=ls, marker=mark, ms=lwg + 1.6, mew=0,
                            solid_capstyle="round", dash_capstyle="round", solid_joinstyle="round",
                            zorder=3 + lw - group)
                    if trials > 1:
                        ax.fill_between(g[0], g[2], g[3], color=c, alpha=alpha, lw=0, zorder=2 + lw - group)
                    other_n |= group == 1
                    ys += list(g[1])

        size_axis(ax, all_sizes)
        ax.set_yscale("log")
        lo, hi = min(ys), max(ys)
        span = (hi / lo) ** 0.04
        ax.set_ylim(lo / span, hi * span)
        ax.yaxis.set_major_locator(FixedLocator(log_ticks(lo / span, hi * span)))
        ax.yaxis.set_minor_locator(NullLocator())
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.3g}"))
        ax.yaxis.set_minor_formatter(NullFormatter())

        sample = (4 if prec == "f32" else 8) * (2 if xform == "cplx" else 1)
        for level, size in caches:
            n = size / (4 * sample)  # working set of one transform: 4 arrays of n samples
            if all_sizes[0] <= n <= all_sizes[-1]:
                cached = True
                right = n > all_sizes[-1] / 8  # near the right edge: label on the left of the line
                ax.axvline(n, color=AXIS, lw=0.9, zorder=1)
                ax.annotate(f"{level} {bytes_label(size)}", (n, 0), xycoords=("data", "axes fraction"),
                            xytext=(-3 if right else 3, 3), textcoords="offset points",
                            ha="right" if right else "left", va="bottom", fontsize=7.5, color=MUTED)
        ax.set_title(f"{prec} {PANEL_TITLES[xform]}", loc="left", fontsize=10.5, fontweight="bold",
                     color=INK, pad=6)
        if i == len(precs) - 1 or (precs[i + 1], xform) not in present:
            ax.set_xlabel("transform length n", fontsize=8.5, color=INK2)
            ax.tick_params(labelbottom=True)
        if j == 0:
            ax.set_ylabel("flops per cycle" if field == "cycles" else "GFLOPS", fontsize=8.5, color=INK2)

    handles = [Line2D([], [], color=color[e], lw=width[e], label=label_for(e)) for e in engines]
    keys = [Line2D([], [], color=MUTED, lw=1.8, label="n a power of two")]
    if other_n:
        keys.append(Line2D([], [], color=MUTED, lw=1.2, ls=(0, (3.5, 2)), label="other n (3·2ᵏ, 5·2ᵏ, ...)"))
    only = [label_for(e) for e in engines if e in POW2_ONLY]
    if only and other_n:
        keys.append(Line2D([], [], color=MUTED, lw=0, label=f"{' and '.join(only)}: powers of two only"))
    if trials > 1:
        keys.append(Patch(color=MUTED, alpha=0.25, lw=0, label=f"range over {trials} trials"))
    if cached:
        keys.append(Line2D([], [], color=AXIS, lw=0.9, label="cache filled by 4 arrays of n samples"))
    if len(files) > 1:
        keys += [Line2D([], [], color=INK2, lw=0, marker=FILE_MARKS[k % len(FILE_MARKS)], ms=5,
                        label=f"{name}: {machine(meta)}") for k, (name, meta, _) in enumerate(files)]
    layout.key(keys, key_ncol)

    title = (f"{label_for(engines[0])} speed" if single else "FFT speed") + "; higher is faster"
    if field == "cycles":
        sub = ("flops per CPU cycle: 5·n·log₂n per complex transform, 2.5·n·log₂n per real one, "
               "÷ cycles per transform")
    else:
        sub = "GFLOPS: 5·n·log₂n per complex transform, 2.5·n·log₂n per real one, ÷ time per transform"
    sub += f"; fastest batch{', median of %d trials' % trials if trials > 1 else ''}"
    if len(files) == 1:
        sub = f"{machine(files[0][1])}  ·  " + sub
    layout.header(title, sub, [] if single else handles, ncol)
    save(args, fig, plt)


def plot_compare(args, files, plt) -> None:
    """B vs A, from a CSV of compare.py."""
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FuncFormatter, MaxNLocator

    if len(files) != 1:
        sys.exit("plot.py: a CSV of compare.py is plotted alone")
    _, meta, rows = files[0]
    a, b = "pffft@A", "pffft@B"
    present, precs, xforms = panels(rows)
    layout = Layout(plt, precs, xforms, 0, 1)
    fig, axes = layout.fig, layout.axes
    for (i, prec), (j, xform) in itertools.product(enumerate(precs), enumerate(xforms)):
        ax = axes[i][j]
        style_axes(ax)
        if (prec, xform) not in present:
            ax.set_visible(False)
            continue
        sizes = sorted({r["n"] for r in rows if (r["prec"], r["xform"]) == (prec, xform)})
        pts = []
        for n in sizes:
            same = [r for r in rows if (r["prec"], r["xform"], r["n"]) == (prec, xform, n)]
            field = "cycles" if all(r["cycles"] is not None for r in same if r["engine"] in (a, b)) else "ns"
            va = [r[field] for r in same if r["engine"] == a and r[field] is not None]
            vb = [r[field] for r in same if r["engine"] == b and r[field] is not None]
            if va and len(va) == len(vb):
                pts.append((n, *results.summary([(y / x - 1) * 100 for x, y in zip(va, vb)])))
        if not pts:
            continue
        ax.axhline(0, color=INK2, lw=1, zorder=2)
        for n, med, lo, hi in pts:
            c = FASTER if hi < 0 else SLOWER if lo > 0 else MUTED
            mark = "o" if results.is_pow2(n) else "^"
            ax.vlines(n, lo, hi, color=c, lw=2, alpha=0.35, zorder=3, capstyle="round")
            ax.plot([n], [med], marker=mark, ms=7 if mark == "^" else 6, color=c, mec=SURFACE, mew=1.2, zorder=4)
        best = min(pts, key=lambda p: p[1])
        worst = max(pts, key=lambda p: p[1])
        for n, med, lo, hi in {best, worst}:
            if abs(med) >= 1:
                ax.annotate(f"{med:+.1f}%".replace("-", "\u2212"), (n, med), xytext=(7, 0), textcoords="offset points",
                            va="center", fontsize=8, color=INK, fontweight="bold")
        size_axis(ax, sizes)
        lo = min(p[2] for p in pts)
        hi = max(p[3] for p in pts)
        pad = 0.08 * max(hi - lo, 2)
        ax.set_ylim(min(lo, 0) - pad, max(hi, 0) + pad)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=7, steps=[1, 2, 2.5, 5, 10]))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: "0" if v == 0 else f"{v:+g}%".replace("-", "\u2212")))
        ax.set_title(f"{prec} {PANEL_TITLES[xform]}", loc="left", fontsize=10.5, fontweight="bold",
                     color=INK, pad=6)
        if i == len(precs) - 1:
            ax.set_xlabel("transform length n", fontsize=8.5, color=INK2)
            ax.tick_params(labelbottom=True)
        if j == 0:
            ax.set_ylabel("change in cost, B vs A", fontsize=8.5, color=INK2)

    keys = [Line2D([], [], color=FASTER, ls="", marker="o", ms=6, label="B faster (whole range below 0)"),
            Line2D([], [], color=SLOWER, ls="", marker="o", ms=6, label="B slower (whole range above 0)"),
            Line2D([], [], color=MUTED, ls="", marker="o", ms=6, label="no clear change"),
            Line2D([], [], color=MUTED, ls="", marker="^", ms=6.5, label="n not a power of two"),
            Line2D([], [], color=MUTED, lw=2, alpha=0.5, label="range over launches and trials")]
    layout.key(keys, 5)
    name = {k: meta.get(k, k) for k in "AB"}
    unit = "CPU cycles" if all(r["cycles"] is not None for r in rows) else "time"
    layout.header(f"pffft {name['B']} vs {name['A']}",
           f"{machine(meta)}  ·  {unit} per pffft_transform_ordered, B relative to A; below 0: B is faster; "
           f"dot: median of {meta.get('launches', '?')} launches × {int(meta.get('trials', 0)) // max(1, int(meta.get('launches', 1)))} trials",
                  [], 1)
    save(args, fig, plt)


def save(args, fig, plt) -> None:
    if args.o:
        fig.savefig(args.o, dpi=200, facecolor=SURFACE)
    else:
        plt.show()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("csv", nargs="+", help="bench_fft CSV files, or one CSV of compare.py")
    ap.add_argument("-o", metavar="FILE", help="write the figure (.png, .svg, ...) instead of showing it")
    ap.add_argument("--engines", metavar="LIST", help="plot only these engines (comma separated)")
    ap.add_argument("--gflops", action="store_true", help="plot GFLOPS (from time) even where cycles were counted")
    args = ap.parse_args()
    files = load(args.csv)
    if args.engines:
        keep = set(args.engines.split(","))
        files = [(name, meta, [r for r in rows if r["engine"] in keep]) for name, meta, rows in files]
        if not any(rows for _, _, rows in files):
            sys.exit(f"plot.py: --engines {args.engines}: none of these engines is in the data")
    try:
        import matplotlib
        if args.o:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        sys.exit("plot.py needs matplotlib 3.7 or later: pip install matplotlib, or run it with\n"
                 "  uv run --with matplotlib python3 bench/plot.py ...")
    plt.rcParams.update({"font.family": "sans-serif",
                         "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
                         "svg.fonttype": "none"})
    if any("A" in meta and "B" in meta for _, meta, _ in files):
        plot_compare(args, files, plt)
    else:
        plot_engines(args, files, plt)


if __name__ == "__main__":
    main()
