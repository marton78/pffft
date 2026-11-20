#!/usr/bin/env python3
"""Plot bench_fft results.

    python3 benchmarks/plot.py results.csv... [--engines LIST] [--time | --per-cycle] [-o out.png|out.svg]

A panel per (precision, transform kind), with n on the x axis:

- for bench_fft CSVs, the speed of each engine, higher is faster: GFLOPS,
  counting 5·n·log2(n) flops per complex transform and 2.5·n·log2(n) per
  real one, as benchFFT does. Where cycles were counted, the time is the
  cycles at the median clock of the run, which keeps the noise of the clock
  out and still gives GFLOPS that compare across devices; elsewhere it is
  the measured time. --time takes the measured time throughout,
  --per-cycle plots flops per CPU cycle. A solid line joins the powers of
  two, a dashed one the other lengths; engines that do only powers of two
  (vdsp, green) get there, dotted, the flops of n at the cost of the next
  power of two, as if their user padded;
- for CSVs that compare.py wrote (-o): B's speed-up over A in percent,
  higher is better, the median and the range over all launches and trials;
  several files (devices) share the panels, one color each. --names A,B
  names the revisions.

Several files (two devices, say) are overlaid, one marker shape per file.
Without -o the figure opens in a window. Needs matplotlib 3.7 or later.
"""
from __future__ import annotations

import argparse
import itertools
import math
import re
import statistics
import sys
from pathlib import Path

import results

POW2_ONLY = {"vdsp", "green"}  # engines that support powers of two only
# An engine keeps its color in every figure. The seven libraries take seven
# colors picked to be as different as possible, in OKLab, for every pair:
# worst pair 17.4 apart (x100) to normal vision, 8.8 under simulated protan or
# deutan vision. fftw-meas is fftw-est at 45% HSV saturation, to read as one
# library.
# pffft@LABEL engines take the colors of the libraries that are absent.
SERIES = {"pffft": "#2a78d6", "fftw-est": "#f06a2e", "fftw-meas": "#f0a584", "vdsp": "#43a26d",
          "pocketfft": "#9a5505", "kissfft": "#e678c0", "fftpack": "#bb1664", "green": "#7a2eb4",
          "mkl": "#898781"}
SPARE = ["#bb1664", "#43a26d", "#7a2eb4", "#9a5505", "#e678c0", "#f06a2e"]
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


DEVICES = {"iPhone12,3": "iPhone 11 Pro (A13)", "iPhone12,1": "iPhone 11 (A13)", "MT6878": "MediaTek MT6878",
           "Cobalt 100": "Microsoft Cobalt 100 (Arm Neoverse N2)"}


def device(cpu):
    """A readable name for a CPU or device id: a known name, or the id
    without trademark signs and clock ("Intel(R) Core(TM) i7-8700B CPU @
    3.20GHz": "Intel Core i7-8700B")."""
    if not cpu or cpu in DEVICES:
        return DEVICES.get(cpu, cpu)
    cpu = re.sub(r"\((R|TM|tm|r)\)", "", cpu)
    cpu = re.sub(r"\s+(CPU\s+)?@\s*[\d.]+\s*GHz", "", cpu)
    return re.sub(r"\s+", " ", cpu).strip()


def system(meta) -> str:
    """The operating system, short: Linux without its kernel version."""
    os_name = meta.get("os") or ""
    return "Linux" if os_name.startswith("Linux") else os_name


def machine(meta) -> str:
    """The device, by a readable name where known, and its system."""
    return ", ".join(x for x in (device(meta.get("cpu")), meta.get("os")) if x) or "unknown machine"


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


def run_clock(rows):
    """The typical clock of a run in GHz: the median of the clock of its
    batches (column ghz) over the rows, or of cycles per ns in files without
    that column; None without cycle counts."""
    clocks = [r["ghz"] for r in rows if r.get("ghz")]
    if not clocks:
        clocks = [r["cycles"] / r["ns"] for r in rows if r["cycles"] and r["ns"]]
    return statistics.median(clocks) if clocks else None


def run_costs(rows, mode):
    """({(engine, prec, xform, n): [cost of each trial]}, clock) for one file.
    mode "clock": ns, from cycles at the run's median clock where the trial
    counted them, else measured (cycles keep the clock's noise out, the
    clock gives GFLOPS that compare across devices); "time": measured ns;
    "cycles": CPU cycles."""
    clock = run_clock(rows) if mode == "clock" else None
    cases = {}
    for r in rows:
        if mode == "cycles":
            cost = r["cycles"]
        elif mode == "clock" and clock and r["cycles"] is not None:
            cost = r["cycles"] / clock
        else:
            cost = r["ns"]
        if cost:
            cases.setdefault((r["engine"], r["prec"], r["xform"], r["n"]), []).append(cost)
    return cases, clock


def panels(every):
    present = {(r["prec"], r["xform"]) for r in every}
    precs = [p for p in ("f32", "f64") if any(p == q for q, _ in present)]
    xforms = [x for x in ("real", "cplx") if any(x == y for _, y in present)]
    return present, precs, xforms


class Layout:
    """A grid of panels under a title, a subtitle and a legend, with a key
    below: sizes in inches, so that text keeps its room at any grid size."""
    PANEL_W, PANEL_H, LEFT, RIGHT, ROW = 6.0, 3.6, 1.05, 0.3, 0.23

    def __init__(self, plt, precs, xforms, legend_rows: int, key_rows: int, right: float = None):
        rows, cols = len(precs), len(xforms)
        self.RIGHT = self.RIGHT if right is None else right
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
    """Speed against n: GFLOPS, or flops per cycle with --per-cycle."""
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter, NullLocator

    every = [r for _, _, rows in files for r in rows]
    mode = "time" if args.time else "cycles" if args.per_cycle else "clock"
    if mode == "cycles" and not any(r["cycles"] is not None for r in every):
        sys.exit("plot.py: --per-cycle: these files have no cycle counts")
    costs = [run_costs(rows, mode) for _, _, rows in files]  # per file: (cases, clock)
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
    labelled = len(files) > 1 and len(engines) == 1  # file names at the line ends need room
    layout = Layout(plt, precs, xforms, 0 if single else math.ceil(len(engines) / ncol), key_rows,
                    right=0.8 if labelled else None)
    fig, axes = layout.fig, layout.axes
    trials = max(len(v) for cases, _ in costs for v in cases.values())
    same = len({tuple(m.get(key) for _, key in CACHES) for _, m, _ in files}) == 1
    caches = [(lv, int(files[0][1][key])) for lv, key in CACHES if same and files[0][1].get(key, "").isdigit()]
    cached = other_n = padded = False
    DASHED, DOTTED = (0, (3.5, 2)), (0, (1, 1.8))

    for (i, prec), (j, xform) in itertools.product(enumerate(precs), enumerate(xforms)):
        ax = axes[i][j]
        style_axes(ax)
        if (prec, xform) not in present:
            ax.set_visible(False)
            continue
        all_sizes = sorted({r["n"] for r in every if (r["prec"], r["xform"]) == (prec, xform)})
        ys = []
        for k, (_, meta, rows) in enumerate(files):
            cases = costs[k][0]
            sizes = sorted({r["n"] for r in rows if (r["prec"], r["xform"]) == (prec, xform)})
            mark = FILE_MARKS[k % len(FILE_MARKS)]
            for e in engines[::-1]:
                pts = ([], [])  # powers of two, other n
                for n in sizes:
                    # an engine for powers of two only pads n up to the next one:
                    # the flops of n at the cost of the padded transform
                    m = n if results.is_pow2(n) or e not in POW2_ONLY else 1 << (n - 1).bit_length()
                    v = cases.get((e, prec, xform, m))
                    if not v:
                        continue
                    # floating-point log2, see commit 419a7e60cd1
                    flops = (5 if xform == "cplx" else 2.5) * n * math.log2(n)
                    vals = [flops / x for x in v]  # per cycle, or per ns: GFLOPS
                    pts[0 if results.is_pow2(n) else 1].append((n, *results.summary(vals)))
                c, lw = color[e], width[e]
                # other n: dashed where the engine does them, dotted where it pads
                other = DOTTED if e in POW2_ONLY else DASHED
                for group, ls, alpha in ((0, "-", 0.14), (1, other, 0.08)):
                    if not pts[group]:
                        continue
                    g = list(zip(*pts[group]))
                    lwg = lw if group == 0 else max(1.1, lw * 0.6)
                    ax.plot(g[0], g[1], color=c, lw=lwg, ls=ls, marker=mark, ms=lwg + 1.6, mew=0,
                            solid_capstyle="round", dash_capstyle="round", solid_joinstyle="round",
                            zorder=3 + lw - group)
                    if trials > 1:
                        ax.fill_between(g[0], g[2], g[3], color=c, alpha=alpha, lw=0, zorder=2 + lw - group)
                    other_n |= group == 1 and e not in POW2_ONLY
                    padded |= group == 1 and e in POW2_ONLY
                    ys += list(g[1])
                if labelled and pts[0]:
                    # files tell apart by marker only: name each at the end of its line
                    n_end, y_end = pts[0][-1][:2]
                    ax.annotate(files[k][0], (n_end, y_end), xytext=(6, 0), textcoords="offset points",
                                va="center", fontsize=8, color=INK2, annotation_clip=False)

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
            ax.set_ylabel("flops per cycle" if mode == "cycles" else "GFLOPS", fontsize=8.5, color=INK2)

    handles = [Line2D([], [], color=color[e], lw=width[e], label=label_for(e)) for e in engines]
    keys = [Line2D([], [], color=MUTED, lw=1.8, label="n a power of two")]
    if other_n:
        keys.append(Line2D([], [], color=MUTED, lw=1.2, ls=DASHED, label="n not a power of two"))
    if padded:
        only = " and ".join(label_for(e) for e in engines if e in POW2_ONLY)
        keys.append(Line2D([], [], color=MUTED, lw=1.2, ls=DOTTED,
                           label=f"n not a power of two, padded to the next one ({only})"))
    if trials > 1:
        keys.append(Patch(color=MUTED, alpha=0.25, lw=0, label=f"range over {trials} trials"))
    if cached:
        keys.append(Line2D([], [], color=AXIS, lw=0.9, label="cache filled by 4 arrays of n samples"))
    if len(files) > 1:
        keys += [Line2D([], [], color=INK2, lw=0, marker=FILE_MARKS[k % len(FILE_MARKS)], ms=5,
                        label=f"{name}: {machine(meta)}") for k, (name, meta, _) in enumerate(files)]
    layout.key(keys, key_ncol)

    title = (f"{label_for(engines[0])} speed" if single else "FFT speed")
    if len(files) == 1:
        title += f" on {device(files[0][1].get('cpu')) or '?'}"
        if system(files[0][1]):
            title += f", {system(files[0][1])}"
    title += "; higher is faster"
    flops = "5·n·log₂n per complex transform, 2.5·n·log₂n per real one"
    clocks = [clock for _, clock in costs]
    if mode == "cycles":
        sub = f"flops per CPU cycle: {flops}, ÷ cycles per transform"
    elif mode == "time" or not any(clocks):
        sub = f"GFLOPS: {flops}, ÷ time per transform"
    elif len(files) == 1:
        sub = f"GFLOPS: {flops}, ÷ CPU cycles per transform at {clocks[0]:.2f} GHz, the median clock of the run"
    else:
        sub = f"GFLOPS: {flops}, ÷ CPU cycles per transform at the median clock of each run where counted, else time"
    sub += f"; fastest batch{', median of %d trials' % trials if trials > 1 else ''}"
    os_name = files[0][1].get("os") or ""
    if len(files) == 1 and os_name != system(files[0][1]):  # the title has it, but not Linux's kernel
        sub = f"{os_name}  ·  " + sub
    layout.header(title, sub, [] if single else handles, ncol)
    save(args, fig, plt)


def plot_compare(args, files, plt) -> None:
    """B vs A from CSVs of compare.py: the change in speed, higher is better,
    one color per file (device)."""
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib.ticker import FuncFormatter, MaxNLocator

    a, b = "pffft@A", "pffft@B"
    names = dict(zip("AB", args.names.split(","))) if args.names else {}
    for k in "AB":
        names.setdefault(k, files[0][1].get(k, k))
    every = [r for _, _, rows in files for r in rows]
    present, precs, xforms = panels(every)
    layout = Layout(plt, precs, xforms, 1, 1)
    fig, axes = layout.fig, layout.axes
    colors = [SERIES["pffft"], "#eb6834", "#1baf7a", "#6250d6"]  # all-pairs separable up to three
    marks = ["o", "s", "D", "^"]
    other_n = False
    for (i, prec), (j, xform) in itertools.product(enumerate(precs), enumerate(xforms)):
        ax = axes[i][j]
        style_axes(ax)
        if (prec, xform) not in present:
            ax.set_visible(False)
            continue
        ax.axhline(0, color=INK2, lw=0.9, zorder=2)
        lows, highs, all_sizes = [0.0], [0.0], set()
        for k, (_, meta, rows) in enumerate(files):
            c, mark = colors[k % len(colors)], marks[k % len(marks)]
            pts = ([], [])
            for n in sorted({r["n"] for r in rows if (r["prec"], r["xform"]) == (prec, xform)}):
                same = [r for r in rows if (r["prec"], r["xform"], r["n"]) == (prec, xform, n)]
                field = "cycles" if all(r["cycles"] is not None for r in same if r["engine"] in (a, b)) else "ns"
                va = [r[field] for r in same if r["engine"] == a and r[field] is not None]
                vb = [r[field] for r in same if r["engine"] == b and r[field] is not None]
                if va and len(va) == len(vb):
                    # speed is flops / cost, so B's speed-up is cost A / cost B
                    pts[0 if results.is_pow2(n) else 1].append(
                        (n, *results.summary([(x / y - 1) * 100 for x, y in zip(va, vb)])))
                    all_sizes.add(n)
            for group, ls, lw, alpha in ((0, "-", 2.0, 0.16), (1, (0, (3.5, 2)), 1.2, 0.09)):
                if not pts[group]:
                    continue
                g = list(zip(*pts[group]))
                ax.plot(g[0], g[1], color=c, lw=lw, ls=ls, marker=mark, ms=4.2 if group == 0 else 3.2,
                        mew=0, solid_capstyle="round", dash_capstyle="round", zorder=4 - group)
                ax.fill_between(g[0], g[2], g[3], color=c, alpha=alpha, lw=0, zorder=3 - group)
                lows.append(min(g[2]))
                highs.append(max(g[3]))
                other_n |= group == 1
        size_axis(ax, sorted(all_sizes))
        lo, hi = min(lows), max(highs)
        pad = 0.06 * max(hi - lo, 2)
        ax.set_ylim(lo - pad, hi + pad)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=7, steps=[1, 2, 2.5, 5, 10]))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: "0" if v == 0 else f"{v:+g}%".replace("-", "\u2212")))
        ax.set_title(f"{prec} {PANEL_TITLES[xform]}", loc="left", fontsize=10.5, fontweight="bold",
                     color=INK, pad=6)
        if i == len(precs) - 1 or (precs[i + 1], xform) not in present:
            ax.set_xlabel("transform length n", fontsize=8.5, color=INK2)
            ax.tick_params(labelbottom=True)
        if j == 0:
            ax.set_ylabel(f"speed of {names['B']} relative to {names['A']}", fontsize=8.5, color=INK2)

    handles = [Line2D([], [], color=colors[k % len(colors)], lw=2, marker=marks[k % len(marks)], ms=4.5, mew=0,
                      label=machine(meta)) for k, (_, meta, _) in enumerate(files)]
    keys = [Line2D([], [], color=MUTED, lw=1.8, label="n a power of two")]
    if other_n:
        keys.append(Line2D([], [], color=MUTED, lw=1.2, ls=(0, (3.5, 2)), label="n not a power of two"))
    launches = files[0][1].get("launches", "?")
    per = int(files[0][1].get("trials", 0)) // max(1, int(files[0][1].get("launches", 1) or 1))
    keys.append(Patch(color=MUTED, alpha=0.25, lw=0, label=f"range over {launches} launches × {per} trials"))
    layout.key(keys, 3)
    unit = "CPU cycles" if all(r["cycles"] is not None for r in every) else "time"
    title = f"pffft {names['B']} vs {names['A']}: speed-up; higher is better"
    sub = (f"change in transforms per {'cycle' if unit == 'CPU cycles' else 'second'} "
           f"(pffft_transform_ordered), {names['B']} relative to {names['A']}, measured side by side in one "
           f"process; median of {launches} launches × {per} trials")
    if len(files) == 1:
        sub = f"{machine(files[0][1])}  ·  " + sub
    layout.header(title, sub, handles if len(files) > 1 else [], max(1, len(files)))
    save(args, fig, plt)


def save(args, fig, plt) -> None:
    if args.o:
        fig.savefig(args.o, dpi=args.dpi, facecolor=SURFACE)
    else:
        plt.show()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("csv", nargs="+", help="bench_fft CSV files, or one CSV of compare.py")
    ap.add_argument("-o", metavar="FILE", help="write the figure (.png, .svg, ...) instead of showing it")
    ap.add_argument("--engines", metavar="LIST", help="plot only these engines (comma separated)")
    ap.add_argument("--dpi", type=int, default=200, help="resolution of a raster -o file (default: 200)")
    ap.add_argument("--names", metavar="A,B", help="names of the two revisions of compare.py CSVs, e.g. 1.1.0,1.3.1")
    units = ap.add_mutually_exclusive_group()
    units.add_argument("--time", action="store_true",
                       help="GFLOPS from the measured time, even where cycles were counted")
    units.add_argument("--per-cycle", action="store_true", help="flops per CPU cycle instead of GFLOPS")
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
                 "  uv run --with matplotlib python3 benchmarks/plot.py ...")
    plt.rcParams.update({"font.family": "sans-serif",
                         "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
                         "svg.fonttype": "none"})
    if any("A" in meta and "B" in meta for _, meta, _ in files):
        plot_compare(args, files, plt)
    else:
        plot_engines(args, files, plt)


if __name__ == "__main__":
    main()
