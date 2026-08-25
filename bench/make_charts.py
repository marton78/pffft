#!/usr/bin/env python3
"""Generate benchmark charts from long-format benchmark samples.

Reads pffft-bench-samples v2 CSVs directly (schema in bench/samples.py);
MFLOPS is derived per repetition and medianed over reps -- never trusted
from the file.

Usage:
    python3 make_charts.py <dir1> [dir2 ...]

Each directory argument may be a colon-separated chain of directories that
are merged into one variant group.  When multiple groups are given, series
labels are suffixed with the directory basename and colors form a gradient
across groups for comparison.

    python3 make_charts.py --evolution CHAIN.json [--target local]

Renders the optimization-chain evolution: (a) per-panel curve charts of the
chain head's measured series and (b) a waterfall of accepted steps' median
Hodges-Lehmann shift_pct per (algo, xform), green when faster / red when
slower.

Charts are written under bench_results/ at the repo root.
"""

import argparse
import json
import os
import sys
from pathlib import Path
from statistics import median

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Patch

_BENCH_DIR = Path(__file__).resolve().parent
if str(_BENCH_DIR) not in sys.path:
    sys.path.insert(0, str(_BENCH_DIR))

from samples import SCHEMA_MAGIC, read_samples  # noqa: E402

REPO_ROOT = _BENCH_DIR.parent
RESULTS_DIR = REPO_ROOT / 'bench_results'


# ---------------------------------------------------------------------------
# Color scheme
# ---------------------------------------------------------------------------

PRODUCT_COLORS = {
    'fftpack':  '#6b7280',   # gray
    'vdsp':     '#e91e8f',   # pink
    'green':    '#16a34a',   # green
    'kiss':     '#ca8a04',   # amber
    'pocket':   '#dc2626',   # red
    'ffts':     '#0891b2',   # teal
    'avfft':    '#059669',   # emerald green
    'fftw':     '#f59e0b',   # orange (default)
    'mkl':      '#9333ea',   # violet
    'pffft':    '#2563eb',   # bold blue
    'pffftu':   '#7c3aed',   # bold purple
}

# FFTW variant-specific orange shades
FFTW_VARIANT_COLORS = {
    'estim': '#f59e0b',
    'auto':  '#ea580c',
    'meas':  '#c2410c',
}

# Display names for products (title-cased)
DISPLAY_NAMES = {
    'fftpack': 'FFTPack',
    'vdsp':    'vDSP',
    'green':   'Green',
    'kiss':    'Kiss',
    'pocket':  'Pocket',
    'ffts':    'FFTS',
    'avfft':   'FFmpeg AVTx',
    'fftw':    'FFTW',
    'mkl':     'MKL',
    'pffft':   'PFFFT',
    'pffftu':  'PFFFT-U',
}

# Draw order: competitors first (thin), then PFFFT on top (thick)
DRAW_ORDER_PRIORITY = {
    'pffft':  100,
    'pffftu': 90,
}


def is_pow2(n):
    """Check if n is a power of two."""
    return n > 0 and (n & (n - 1)) == 0


def pow2_mask(sizes):
    """Return a boolean numpy array indicating which sizes are powers of two."""
    arr = np.array(sizes, dtype=np.int64)
    return (arr > 0) & ((arr & (arr - 1)) == 0)


def split_algo(algo):
    """Split an algo id like 'fftw-estim' into (product, variant)."""
    product, _, variant = algo.partition('-')
    return product, (variant or 'default')


def read_provenance(path):
    """Read just the provenance header of a samples file.

    Returns dict, or None when the file lacks the samples-v2 magic.
    """
    try:
        with open(path) as f:
            if f.readline().rstrip('\n') != SCHEMA_MAGIC:
                return None
            prov = {}
            for line in f:
                line = line.rstrip('\n')
                if not line.startswith('#'):
                    break
                k, sep, v = line[2:].partition('=')
                if sep:
                    prov[k] = v.strip()
            return prov
    except OSError:
        return None


def read_samples_file(path):
    """Read one long-format samples file.

    Returns {product: {(prec, xform): (sizes[], median_mflops[])}}, where
    product is the row algo id ('pffft', 'fftw-estim', ...) and each value
    medians derived MFLOPS over repetitions.  Returns None when path is not
    a readable pffft-bench-samples v2 file.
    """
    try:
        _, rows = read_samples(path)
    except (OSError, ValueError):
        return None
    acc = {}  # algo -> (prec, xform) -> size -> [mflops]
    for s in rows:
        m = s.mflops
        if m > 0:
            acc.setdefault(s.algo, {}).setdefault(
                (s.prec, s.xform), {}).setdefault(s.size, []).append(m)
    out = {}
    for algo, panels in acc.items():
        entry = {}
        for key, by_size in panels.items():
            sizes = sorted(by_size)
            entry[key] = (sizes, [median(by_size[sz]) for sz in sizes])
        out[algo] = entry
    return out


def get_color(product, variant, dir_index=None, num_dirs=1):
    """Return color for a product-variant combination.

    When num_dirs > 1, generates a gradient for each product across
    directories so that successive versions are visually distinct.
    """
    if num_dirs > 1 and dir_index is not None:
        return _gradient_color(product, dir_index, num_dirs)
    if product == 'fftw' and variant in FFTW_VARIANT_COLORS:
        return FFTW_VARIANT_COLORS[variant]
    return PRODUCT_COLORS.get(product, '#888888')


# Per-product gradient endpoints: (start_rgb, end_rgb)
# Light/cool -> dark/warm so early = faint, latest = bold
_PRODUCT_GRADIENTS = {
    'pffft':  ((0.68, 0.85, 0.96), (0.08, 0.20, 0.70)),   # light sky -> dark blue
    'pffftu': ((0.82, 0.70, 0.96), (0.45, 0.08, 0.70)),   # light lavender -> deep purple
    'fftpack': ((0.85, 0.85, 0.85), (0.35, 0.38, 0.42)),
    'fftw':    ((0.99, 0.82, 0.55), (0.76, 0.25, 0.05)),
    'vdsp':    ((0.96, 0.70, 0.85), (0.76, 0.08, 0.40)),
    'green':   ((0.70, 0.92, 0.70), (0.05, 0.50, 0.18)),
    'kiss':    ((0.95, 0.88, 0.60), (0.60, 0.42, 0.02)),
    'pocket':  ((0.95, 0.70, 0.70), (0.70, 0.10, 0.10)),
    'avfft':   ((0.60, 0.95, 0.80), (0.01, 0.50, 0.35)),   # light mint -> deep emerald
    'mkl':     ((0.82, 0.70, 0.96), (0.50, 0.15, 0.80)),
}


def _gradient_color(product, dir_index, num_dirs):
    """Interpolate between gradient endpoints for a product."""
    grad = _PRODUCT_GRADIENTS.get(product)
    if grad is None:
        t = dir_index / max(num_dirs - 1, 1)
        grey = 0.75 - 0.4 * t
        return (grey, grey, grey)
    start, end = grad
    t = dir_index / max(num_dirs - 1, 1)
    r = start[0] + t * (end[0] - start[0])
    g = start[1] + t * (end[1] - start[1])
    b = start[2] + t * (end[2] - start[2])
    return (r, g, b)


def make_label(product, variant, dir_suffix=None):
    """Build a human-readable series label."""
    display = DISPLAY_NAMES.get(product, product)
    if variant and variant != 'default':
        label = f'{display} {variant}'
    else:
        label = display
    if dir_suffix:
        label = f'{label} ({dir_suffix})'
    return label


def draw_order_key(product):
    """Sort key: higher priority products are drawn last (on top)."""
    return DRAW_ORDER_PRIORITY.get(product, 0)


def plot_series(ax, sizes, mflops, label, color, linewidth=1.5,
                markersize=3, alpha=0.8, linestyle='-'):
    """Plot one series, filtering to sizes >= 32 with positive mflops.

    Marker shape is determined per point: circle for pow2, triangle for non-pow2.
    """
    xs, ys = [], []
    for s, m in zip(sizes, mflops):
        if m > 0 and s >= 32:
            xs.append(s)
            ys.append(m)
    if not xs:
        return
    xs_arr = np.array(xs)
    ys_arr = np.array(ys)
    mask = pow2_mask(xs)
    # Draw the connecting line without markers, carrying the legend label
    ax.plot(xs_arr, ys_arr, color=color, linewidth=linewidth, alpha=alpha,
            linestyle=linestyle, label=label)
    # Overlay pow2 points as circles
    if mask.any():
        ax.plot(xs_arr[mask], ys_arr[mask], 'o', color=color,
                markersize=markersize, alpha=alpha)
    # Overlay non-pow2 points as triangles
    if (~mask).any():
        ax.plot(xs_arr[~mask], ys_arr[~mask], 'x', color=color,
                markersize=markersize, alpha=alpha)


def make_chart(ax, series_list, title, num_variants=1):
    """Draw a single chart panel.

    series_list: [(product, variant, label, sizes, mflops, var_index), ...]
    """
    # Sort so competitors are drawn first, PFFFT on top
    ordered = sorted(series_list, key=lambda s: draw_order_key(s[0]))

    has_np2 = False
    for product, variant, label, sizes, mflops, var_index in ordered:
        color = get_color(product, variant, var_index, num_variants)
        is_pffft = product in ('pffft', 'pffftu')
        lw = 2.2 if is_pffft else 1.2
        ms = 6 if is_pffft else 5
        al = 1.0 if is_pffft else 0.85
        plot_series(ax, sizes, mflops, label, color,
                    linewidth=lw, markersize=ms, alpha=al)
        if not has_np2:
            mask = pow2_mask([s for s, m in zip(sizes, mflops) if m > 0 and s >= 32])
            has_np2 = (~mask).any()

    ax.set_xscale('log', base=2)
    ax.set_yscale('log', base=10)
    ax.set_xlabel('FFT size', fontsize=11)
    ax.set_ylabel('MFlops (higher = better)', fontsize=11)
    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, which='both')
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(
        lambda x, _: f'{int(x):,}' if x >= 1 else ''))

    # Add marker-type legend entries when non-pow2 data is present
    handles, labels = ax.get_legend_handles_labels()
    if has_np2:
        from matplotlib.lines import Line2D
        handles += [
            Line2D([0], [0], marker='o', color='gray', linestyle='none',
                   markersize=5, label='pow2'),
            Line2D([0], [0], marker='x', color='gray', linestyle='none',
                   markersize=5, label='non-pow2'),
        ]
    ax.legend(handles=handles, fontsize=8, loc='best', framealpha=0.9)


PRECISION_LABELS = {
    'flt': 'Single-precision',
    'dbl': 'Double-precision',
}
TRANSFORM_LABELS = {
    'real': 'real FFT',
    'cplx': 'complex FFT',
}


def panel_title(precision, transform, prov=None):
    """Build a chart panel title, including dominant provenance."""
    prec = PRECISION_LABELS.get(precision, precision)
    xform = TRANSFORM_LABELS.get(transform, transform)
    title = f'{prec} {xform}'
    if prov:
        tag = ' @ '.join(x for x in (prov.get('label'), prov.get('host'))
                         if x)
        if tag:
            title = f'{title} \u2014 {tag}'
    return title


def scan_directory(dirpath):
    """Scan a directory for long-format samples CSV files.

    Globs *.csv; files without the samples-v2 magic are skipped gracefully.
    Products (algos) are merged across all readable files in the directory.

    Returns (provenance, panels): provenance is the header dict of the first
    readable file (the directory's dominant provenance) or None; panels maps
    (prec, xform) -> [(product, variant, sizes, mflops)].
    """
    panels = {}
    dominant_prov = None
    try:
        entries = sorted(os.listdir(dirpath))
    except OSError:
        return dominant_prov, panels

    for fname in entries:
        if not fname.endswith('.csv'):
            continue
        path = os.path.join(dirpath, fname)
        prov = read_provenance(path)
        if prov is None:
            print(f'make_charts: skipping {path}: not a '
                  f'{SCHEMA_MAGIC!r} file', file=sys.stderr)
            continue
        if dominant_prov is None:
            dominant_prov = prov
        data = read_samples_file(path)
        for algo, algo_panels in (data or {}).items():
            product, variant = split_algo(algo)
            for key, (sizes, mflops) in algo_panels.items():
                panels.setdefault(key, []).append(
                    (product, variant, sizes, mflops))

    return dominant_prov, panels


def merge_dirs(dirpaths):
    """Merge scan results from multiple colon-chained directories.

    For each (prec, xform, product, variant), takes the per-size median over
    every contributing file's medians.  Returns (provenance, panels) with
    provenance taken from the first directory that yields one.
    """
    merged = {}  # (prec, xform) -> {(product, variant) -> size -> [mflops]}
    prov = None
    for dirpath in dirpaths:
        dprov, dpanels = scan_directory(dirpath)
        if prov is None:
            prov = dprov
        for key, series in dpanels.items():
            for product, variant, sizes, mflops in series:
                by_size = merged.setdefault(key, {}).setdefault(
                    (product, variant), {})
                for sz, m in zip(sizes, mflops):
                    by_size.setdefault(sz, []).append(m)

    result = {}
    for key, pv_dict in merged.items():
        result[key] = []
        for (product, variant), by_size in pv_dict.items():
            sizes = sorted(by_size)
            result[key].append((
                product, variant,
                sizes,
                [median(by_size[sz]) for sz in sizes],
            ))
    return prov, result


def render_panels(all_panels, prov_list, output_dir, suptitle_base):
    """Write per-panel webp charts plus a combined grid under output_dir."""
    panel_order = [
        ('flt', 'real'),
        ('flt', 'cplx'),
        ('dbl', 'real'),
        ('dbl', 'cplx'),
    ]
    active_panels = [(k, all_panels[k]) for k in panel_order if k in all_panels]
    panel_names = {
        ('flt', 'real'):  'float_real',
        ('flt', 'cplx'):  'float_cplx',
        ('dbl', 'real'):  'double_real',
        ('dbl', 'cplx'):  'double_cplx',
    }
    num_variants = len(prov_list)
    saved = []
    for (prec, xform), series_list in active_panels:
        title = panel_title(prec, xform, next(
            (p for p in prov_list if p), None))
        fig, ax = plt.subplots(figsize=(11, 6.5))
        make_chart(ax, series_list, title, num_variants)
        fig.tight_layout()
        name = panel_names.get((prec, xform), f'{prec}_{xform}')
        outpath = os.path.join(output_dir, f'bench_{name}.webp')
        fig.savefig(outpath, dpi=150, format='webp')
        plt.close(fig)
        saved.append(outpath)
        print(f'Saved {outpath}')
    n = len(active_panels)
    if n >= 2:
        prov = next((p for p in prov_list if p), None)
        suptitle = suptitle_base
        if prov:
            tag = ' @ '.join(x for x in (prov.get('label'), prov.get('host'))
                             if x)
            if tag:
                suptitle = f'{suptitle_base} \u2014 {tag}'
        fig = draw_combined(active_panels, prov, suptitle, num_variants)
        outpath = os.path.join(output_dir, 'bench_all.webp')
        fig.savefig(outpath, dpi=150, format='webp')
        plt.close(fig)
        saved.append(outpath)
        print(f'Saved {outpath}')
    return saved


def draw_combined(active_panels, prov, suptitle, num_variants=1):
    """Draw active panels on a single figure (2x2 grid, smaller when n < 4)."""
    n = len(active_panels)
    rows = 2 if n > 2 else 1
    cols = 2 if n > 1 else 1
    figsize = (11, 6.5) if n == 1 else (20, 12 if rows == 2 else 7)
    fig, axes = plt.subplots(rows, cols, figsize=figsize)
    if not hasattr(axes, '__iter__'):
        axes = [axes]
    flat = [ax for row in axes
            for ax in (row if hasattr(row, '__iter__') else [row])]

    for i, ((prec, xform), series_list) in enumerate(active_panels):
        if i < len(flat):
            make_chart(flat[i], series_list,
                       panel_title(prec, xform, prov), num_variants)

    for j in range(n, len(flat)):
        flat[j].set_visible(False)

    fig.suptitle(suptitle, fontsize=16, fontweight='bold', y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    return fig


def run_directory_mode(dir_args):
    """Original multi-directory comparison mode."""
    # Each argument may be colon-separated paths forming one variant group
    groups = [arg.split(':') for arg in dir_args]
    groups = [[os.path.abspath(d) for d in g] for g in groups]
    multi = len(groups) > 1

    # all_panels: (prec, xform) -> [(product, variant, label, sizes, mflops, var_index)]
    all_panels = {}
    prov_list = []

    for var_index, dirpaths in enumerate(groups):
        prov, panels = merge_dirs(dirpaths)
        prov_list.append(prov)
        dir_suffix = '+'.join(os.path.basename(d) for d in dirpaths) if multi else None

        for key, series in panels.items():
            for product, variant, sizes, mflops in series:
                label = make_label(product, variant, dir_suffix)
                all_panels.setdefault(key, []).append(
                    (product, variant, label, sizes, mflops, var_index))

    if not all_panels:
        print('No benchmark samples found in the given directories!',
              file=sys.stderr)
        sys.exit(1)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    render_panels(all_panels, prov_list, RESULTS_DIR, 'PFFFT Benchmark')


# ---------------------------------------------------------------------------
# Evolution mode
# ---------------------------------------------------------------------------

def resolve_sample_file(p):
    """Resolve a chain sample_files entry to an existing Path, or None."""
    path = Path(p)
    if path.is_absolute():
        return path if path.exists() else None
    for base in (REPO_ROOT, REPO_ROOT / '.perf'):
        cand = base / path
        if cand.exists():
            return cand
    return None


def filter_head_files(files, tree, label):
    """Keep only sample files belonging to the head build.

    `accept` stores base+variant files together in one step's sample_files,
    so the head curve must not aggregate baseline builds.  A file matches
    when its provenance git_tree equals `tree`; when either side lacks a
    git_tree, fall back to a provenance label == `label` match.  Files with
    an unreadable header are dropped (aggregate_samples skips them anyway).
    """
    kept = []
    for p in files:
        prov = read_provenance(p)
        if prov is None:
            continue
        if tree and prov.get('git_tree'):
            if prov['git_tree'] == tree:
                kept.append(p)
        elif label and prov.get('label') == label:
            kept.append(p)
    return kept


def aggregate_samples(files, target=None):
    """Median derived MFLOPS per (prec, xform, algo, size) across files/reps.

    Returns panels: (prec, xform) -> [(product, variant, sizes, mflops)].
    Files whose provenance target differs from `target` are skipped.
    """
    acc = {}  # (prec, xform, algo, size) -> [mflops]
    for p in files:
        try:
            prov, rows = read_samples(p)
        except (OSError, ValueError):
            continue
        if target and prov.get('target') != target:
            continue
        for s in rows:
            m = s.mflops
            if m > 0:
                acc.setdefault((s.prec, s.xform, s.algo, s.size), []).append(m)

    algos = {}  # (prec, xform) -> algo -> size -> mflops
    for (prec, xform, algo, size), ms in acc.items():
        algos.setdefault((prec, xform), {}).setdefault(algo, {})[size] = median(ms)

    panels = {}
    for (prec, xform), by_algo in algos.items():
        lst = []
        for algo, by_size in by_algo.items():
            product, variant = split_algo(algo)
            sizes = sorted(by_size)
            lst.append((product, variant, sizes,
                        [by_size[sz] for sz in sizes]))
        panels[(prec, xform)] = lst
    return panels


def load_evolution_chain(chain_path):
    try:
        with open(chain_path) as f:
            chain = json.load(f)
    except (OSError, ValueError) as e:
        raise SystemExit(f'--evolution: cannot read {chain_path}: {e}')
    if not isinstance(chain, dict):
        raise SystemExit(f'--evolution: {chain_path} is not a chain manifest')
    return chain


def run_evolution_mode(chain_path, target):
    """Render (a) chain-head curves and (b) accepted-steps waterfall."""
    chain = load_evolution_chain(chain_path)
    steps = [s for s in chain.get('steps') or [] if s.get('accepted')]
    base = chain.get('base') or {}

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    saved = []

    # --- (a) chain head series -------------------------------------------
    if steps:
        head_name = steps[-1].get('label') or 'head'
        head_tree = (steps[-1].get('tree') or '')[:8]
        head_full_tree = steps[-1].get('tree') or ''
        raw_files = steps[-1].get('sample_files') or []
    else:
        head_name = base.get('label') or 'base'
        head_tree = (base.get('tree') or '')[:8]
        head_full_tree = base.get('tree') or ''
        raw_files = []

    head_files = []
    for p in raw_files:
        resolved = resolve_sample_file(p)
        if resolved is None:
            print(f'evolution: missing sample file {p}', file=sys.stderr)
        else:
            head_files.append(resolved)
    # accept stores base+variant files together in one step; keep only the
    # files that actually belong to the head build.
    kept = filter_head_files(head_files, head_full_tree, head_name)
    if len(kept) < len(head_files):
        print(f'evolution: dropped {len(head_files) - len(kept)} sample '
              f'file(s) not belonging to head {head_name!r}', file=sys.stderr)
    head_files = kept
    panels = aggregate_samples(head_files, target=target)
    if not panels:
        raise SystemExit(f'evolution: no readable sample data for target '
                         f'{target!r} in chain head {head_name!r}')

    suffix = f' @ {head_tree}' if head_tree else ''
    prov = {'label': f'{head_name}{suffix}'}
    head_panels = []
    for key, series in sorted(panels.items()):
        entries = [(product, variant, make_label(product, variant),
                    sizes, mflops, 0)
                   for product, variant, sizes, mflops in series]
        head_panels.append((key, entries))
    fig = draw_combined(head_panels, prov,
                        f'Evolution head: {head_name}{suffix}')
    outpath = os.path.join(RESULTS_DIR, f'evolution-{target}-curves.png')
    fig.savefig(outpath, dpi=150, format='png')
    plt.close(fig)
    saved.append(outpath)
    print(f'Saved {outpath}')

    # --- (b) waterfall of accepted steps ---------------------------------
    groups = {}  # (algo, prec, xform) -> [(step_label, shift_pct)]
    for step in steps:
        for v in step.get('verdicts') or []:
            shift = v.get('shift_pct')
            if shift is None:
                continue
            if (v.get('target') or 'local') != target:
                continue
            key = (v.get('algo'), v.get('prec'), v.get('xform'))
            groups.setdefault(key, []).append(
                (step.get('label'), float(shift)))

    outpath = os.path.join(RESULTS_DIR, f'evolution-{target}-waterfall.png')
    fig, ax = plt.subplots(figsize=(max(9, 1.15 * len(groups)), 6))
    if groups:
        keys = sorted(groups)
        max_steps = max(len(v) for v in groups.values())
        width = 0.8 / max_steps
        for gi, key in enumerate(keys):
            vals = groups[key]
            for si, (_, pct) in enumerate(vals):
                offset = (si - (len(vals) - 1) / 2) * width
                color = '#16a34a' if pct < 0 else '#dc2626'
                ax.bar(gi + offset, pct, width=width * 0.92, color=color,
                       edgecolor='none')
        ax.set_xticks(range(len(keys)))
        ax.set_xticklabels([f'{a}\n{p}/{x}' for a, p, x in keys], fontsize=8)
        step_labels = [s.get('label') or '?' for s in steps]
        ax.set_title(
            f'Accepted-step HL shift \u2014 target {target} \u2014 steps: '
            + ' \u2192 '.join(step_labels),
            fontsize=12, fontweight='bold')
    else:
        ax.text(0.5, 0.5, f'no accepted steps with verdicts for target '
                          f'{target!r}', ha='center', va='center',
                transform=ax.transAxes)
    ax.axhline(0, color='black', linewidth=0.8)
    ax.set_ylabel('HL shift % (negative = faster)', fontsize=11)
    ax.legend(handles=[
        Patch(facecolor='#16a34a', label='faster (< 0)'),
        Patch(facecolor='#dc2626', label='slower (> 0)'),
    ], fontsize=8, loc='best', framealpha=0.9)
    ax.grid(True, axis='y', alpha=0.3)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150, format='png')
    plt.close(fig)
    saved.append(outpath)
    print(f'Saved {outpath}')
    return saved


def main(argv=None):
    ap = argparse.ArgumentParser(
        description='Render PFFFT benchmark charts from long-format samples.')
    ap.add_argument('dirs', nargs='*',
                    help='result directories (colon chains allowed)')
    ap.add_argument('--evolution', metavar='CHAIN',
                    help='render evolution charts for a bench_chain.json')
    ap.add_argument('--target', default='local',
                    help='target to filter evolution data by (default local)')
    args = ap.parse_args(argv)

    if args.evolution:
        chain_path = args.evolution
        if not os.path.isabs(chain_path):
            cand = Path.cwd() / chain_path
            chain_path = str(cand) if cand.exists() else os.path.abspath(
                os.path.join(str(REPO_ROOT), chain_path))
        run_evolution_mode(chain_path, args.target)
        return

    if not args.dirs:
        ap.error('give result directories or --evolution CHAIN')
    run_directory_mode(args.dirs)


if __name__ == '__main__':
    main()
