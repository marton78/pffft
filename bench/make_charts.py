#!/usr/bin/env python3
"""Generate benchmark charts from per-library CSV results.

Reads per-library CSV files with the naming pattern:
    <product>-<variant>-<flt|dbl>-<real|cplx>.csv

Each CSV has columns: size,prep_ms,num_iter,mflops,duration_sec

Usage:
    python3 make_charts.py <dir1> [dir2 ...]

When multiple directories are given, series labels are suffixed with the
directory basename for comparison.  Output .webp files are written to the
first directory.
"""

import csv
import json
import os
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


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


def parse_csv_filename(fname):
    """Parse a CSV filename into (product_variant, precision, transform).

    Pattern: <product-variant>-<flt|dbl>-<real|cplx>.csv
    Parse from the right: last two segments before .csv are precision and
    transform; everything before that is product-variant.
    """
    base = fname.removesuffix('.csv') if fname.endswith('.csv') else None
    if base is None:
        return None
    parts = base.split('-')
    if len(parts) < 3:
        return None
    transform = parts[-1]
    precision = parts[-2]
    if precision not in ('flt', 'dbl') or transform not in ('real', 'cplx'):
        return None
    product_variant = '-'.join(parts[:-2])
    return product_variant, precision, transform


def read_per_library_csv(path):
    """Read a per-library CSV file.  Returns list of (size, mflops) tuples."""
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                size = int(row['size'])
                mflops = float(row['mflops'])
                rows.append((size, mflops))
            except (ValueError, KeyError):
                continue
    return rows


def get_color(product, variant):
    """Return color for a product-variant combination."""
    if product == 'fftw' and variant in FFTW_VARIANT_COLORS:
        return FFTW_VARIANT_COLORS[variant]
    return PRODUCT_COLORS.get(product, '#888888')


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
                marker='o', markersize=3, alpha=0.8, linestyle='-'):
    """Plot one series, filtering to sizes >= 32 with positive mflops."""
    xs, ys = [], []
    for s, m in zip(sizes, mflops):
        if m > 0 and s >= 32:
            xs.append(s)
            ys.append(m)
    if not xs:
        return
    ax.plot(xs, ys, marker=marker, color=color, label=label,
            markersize=markersize, linewidth=linewidth, alpha=alpha,
            linestyle=linestyle)


def make_chart(ax, series_list, title):
    """Draw a single chart panel.

    series_list: [(product, variant, label, sizes, mflops), ...]
    """
    # Sort so competitors are drawn first, PFFFT on top
    ordered = sorted(series_list, key=lambda s: draw_order_key(s[0]))

    for product, variant, label, sizes, mflops in ordered:
        color = get_color(product, variant)
        is_pffft = product in ('pffft', 'pffftu')
        lw = 2.2 if is_pffft else 1.2
        ms = 4 if is_pffft else 3
        al = 1.0 if is_pffft else 0.85
        plot_series(ax, sizes, mflops, label, color,
                    linewidth=lw, markersize=ms, alpha=al)

    ax.set_xscale('log', base=2)
    ax.set_yscale('log', base=10)
    ax.set_xlabel('FFT size', fontsize=11)
    ax.set_ylabel('MFlops (higher = better)', fontsize=11)
    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, which='both')
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(
        lambda x, _: f'{int(x):,}' if x >= 1 else ''))
    ax.legend(fontsize=8, loc='best', framealpha=0.9)


def load_info(dirpath):
    """Load info.json from a directory, return dict or empty dict."""
    info_path = os.path.join(dirpath, 'info.json')
    if os.path.exists(info_path):
        try:
            with open(info_path) as f:
                return json.load(f)
        except (json.JSONDecodeError, OSError):
            pass
    return {}


PRECISION_LABELS = {
    'flt': 'Single-precision',
    'dbl': 'Double-precision',
}
TRANSFORM_LABELS = {
    'real': 'real FFT',
    'cplx': 'complex FFT',
}


def panel_title(precision, transform, info_list):
    """Build a chart panel title."""
    prec = PRECISION_LABELS.get(precision, precision)
    xform = TRANSFORM_LABELS.get(transform, transform)
    title = f'{prec} {xform}'
    # Add CPU info from first directory that has it
    for info in info_list:
        cpu = info.get('cpu')
        if cpu:
            title = f'{title} \u2014 {cpu}'
            break
    return title


def scan_directory(dirpath):
    """Scan a directory for per-library CSV files.

    Returns dict: (precision, transform) -> [(product, variant, sizes, mflops)]
    """
    panels = {}
    try:
        entries = os.listdir(dirpath)
    except OSError:
        return panels

    for fname in sorted(entries):
        if not fname.endswith('.csv'):
            continue
        parsed = parse_csv_filename(fname)
        if parsed is None:
            continue
        product_variant, precision, transform = parsed
        # Split product from variant: product is the first segment
        pv_parts = product_variant.split('-', 1)
        product = pv_parts[0]
        variant = pv_parts[1] if len(pv_parts) > 1 else 'default'

        data = read_per_library_csv(os.path.join(dirpath, fname))
        if not data:
            continue
        sizes = [d[0] for d in data]
        mflops = [d[1] for d in data]

        key = (precision, transform)
        panels.setdefault(key, []).append((product, variant, sizes, mflops))

    return panels


def main():
    if len(sys.argv) < 2:
        print('Usage: make_charts.py <dir1> [dir2 ...]', file=sys.stderr)
        sys.exit(1)

    dirs = sys.argv[1:]
    multi = len(dirs) > 1
    output_dir = dirs[0]

    # Collect data from all directories
    # all_panels: (prec, xform) -> [(product, variant, label, sizes, mflops)]
    all_panels = {}
    info_list = []

    for dirpath in dirs:
        dirpath = os.path.abspath(dirpath)
        info = load_info(dirpath)
        info_list.append(info)
        dir_suffix = os.path.basename(dirpath) if multi else None

        panels = scan_directory(dirpath)
        for key, series in panels.items():
            for product, variant, sizes, mflops in series:
                label = make_label(product, variant, dir_suffix)
                all_panels.setdefault(key, []).append(
                    (product, variant, label, sizes, mflops))

    if not all_panels:
        print('No CSV data found in the given directories!', file=sys.stderr)
        sys.exit(1)

    # Canonical panel order
    panel_order = [
        ('flt', 'real'),
        ('flt', 'cplx'),
        ('dbl', 'real'),
        ('dbl', 'cplx'),
    ]

    # Keep only panels that have data, in canonical order
    active_panels = [(k, all_panels[k]) for k in panel_order if k in all_panels]

    # --- Individual charts ---
    panel_names = {
        ('flt', 'real'):  'float_real',
        ('flt', 'cplx'):  'float_cplx',
        ('dbl', 'real'):  'double_real',
        ('dbl', 'cplx'):  'double_cplx',
    }
    for (prec, xform), series_list in active_panels:
        title = panel_title(prec, xform, info_list)
        fig, ax = plt.subplots(figsize=(11, 6.5))
        make_chart(ax, series_list, title)
        fig.tight_layout()
        name = panel_names.get((prec, xform), f'{prec}_{xform}')
        outpath = os.path.join(output_dir, f'bench_{name}.webp')
        fig.savefig(outpath, dpi=150, format='webp')
        plt.close(fig)
        print(f'Saved {outpath}')

    # --- Combined 2x2 chart ---
    n = len(active_panels)
    if n >= 2:
        rows = 2 if n > 2 else 1
        cols = 2
        fig, axes = plt.subplots(rows, cols, figsize=(20, 12 if rows == 2 else 7))
        if rows == 1:
            axes = [axes]
        flat = [ax for row in axes
                for ax in (row if hasattr(row, '__iter__') else [row])]

        for i, ((prec, xform), series_list) in enumerate(active_panels):
            if i < len(flat):
                title = panel_title(prec, xform, info_list)
                make_chart(flat[i], series_list, title)

        for j in range(n, len(flat)):
            flat[j].set_visible(False)

        # Suptitle from info
        suptitle = 'PFFFT Benchmark'
        for info in info_list:
            cpu = info.get('cpu')
            if cpu:
                suptitle = f'PFFFT Benchmark \u2014 {cpu}'
                break

        fig.suptitle(suptitle, fontsize=16, fontweight='bold', y=0.98)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        outpath = os.path.join(output_dir, 'bench_all.webp')
        fig.savefig(outpath, dpi=150, format='webp')
        plt.close(fig)
        print(f'Saved {outpath}')


if __name__ == '__main__':
    main()
