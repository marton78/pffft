#!/usr/bin/env python3
"""Plot MFLOPS evolution across library versions from pffft-bench-samples files.

Each input file represents one "version" (its `label` provenance key, or the
filename stem if no label is present). For each requested FFT size, renders
one panel with versions along the x-axis (in the given/discovered order) and
a box or violin plot of that version's MFLOPS distribution across reps.

Usage
-----
    ./bench/plot_evolution.py .perf/samples/*.csv --sizes 256,1024,4096

    ./bench/plot_evolution.py /tmp/evosamples/*.csv \\
        --order master,float-fma,double-rewrite,vmsub-cplxmul,final \\
        --sizes 256,1024,4096 --algo pffft --xform real --prec flt \\
        --kind violin --out bench_results/evolution_pffft_real_flt.png

Reusable as a library too:
    from plot_evolution import load_versions, plot_evolution
"""
from __future__ import annotations

import argparse
import glob
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from samples import read_samples, group_key  # noqa: E402


def load_versions(paths: list[str]) -> dict[str, list]:
    """Read each samples file into {version_label: [Sample, ...]}.

    Later files with the same label extend the earlier ones (useful for
    incrementally-collected data), preserving each label's first-seen order.
    """
    versions: dict[str, list] = {}
    for pattern in paths:
        for path in sorted(glob.glob(pattern)) or [pattern]:
            if not Path(path).is_file():
                continue
            prov, rows = read_samples(path)
            label = prov.get("label") or Path(path).stem
            versions.setdefault(label, []).extend(rows)
    return versions


def select_mflops(rows: list, algo: str, prec: str, xform: str, size: int) -> list[float]:
    return [s.mflops for s in rows
            if s.algo == algo and s.prec == prec and s.xform == xform and s.size == size]


def plot_evolution(versions: dict[str, list], sizes: list[int], *,
                    algo: str = "pffft", prec: str = "flt", xform: str = "real",
                    order: list[str] | None = None, kind: str = "box",
                    title: str | None = None,
                    noop: "set[str] | dict[int, set[str]] | None" = None):
    """Build the figure (import matplotlib lazily so this module stays
    importable without a plotting backend for non-plotting callers).

    `noop` names version labels expected to show NO real effect relative to
    the immediately preceding label on the x-axis (e.g. a double-only change
    panel doesn't exercise). Such steps get a hatched box and a "=" marker
    over the gap from the previous step, so a reader can tell "flat here"
    was predicted, not just observed.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = order or list(versions.keys())
    missing = [lab for lab in labels if lab not in versions]
    if missing:
        raise ValueError(f"--order names version(s) not found in the input data: {missing}")
    noop = noop or {}

    fig, axes = plt.subplots(1, len(sizes), figsize=(5 * len(sizes), 5), squeeze=False)
    axes = axes[0]
    color = "#3b7dd8"
    noop_color = "#b0b0b0"
    any_noop_rendered = False


    for ax, size in zip(axes, sizes):
        data, present = [], []
        for lab in labels:
            vals = select_mflops(versions[lab], algo, prec, xform, size)
            if vals:
                data.append(vals)
                present.append(lab)
        if not data:
            ax.set_title(f"N={size} (no data)")
            ax.axis("off")
            continue

        positions = list(range(1, len(data) + 1))
        noop_here = noop.get(size, set()) if isinstance(noop, dict) else noop
        is_noop = [lab in noop_here for lab in present]
        if kind == "violin":
            parts = ax.violinplot(data, positions=positions, showmedians=True, widths=0.8)
            for body, noopflag in zip(parts["bodies"], is_noop):
                body.set_facecolor(noop_color if noopflag else color)
                body.set_alpha(0.6)
        else:
            bp = ax.boxplot(data, positions=positions, widths=0.6, showfliers=True,
                            patch_artist=True,
                            boxprops=dict(facecolor=color, alpha=0.6),
                            medianprops=dict(color="black", linewidth=1.5))
            for box, noopflag in zip(bp["boxes"], is_noop):
                if noopflag:
                    box.set_facecolor(noop_color)
                    box.set_hatch("//")
                    box.set_edgecolor("#606060")

        # "=" markers over the gap from each no-op step's predecessor, plus
        # a legend proxy so the meaning is explained once per figure.
        ymax = max(max(v) for v in data)
        ymin = min(min(v) for v in data)
        pad = (ymax - ymin) * 0.08 or ymax * 0.05
        for i, noopflag in enumerate(is_noop):
            if noopflag and i > 0:
                mid = (positions[i - 1] + positions[i]) / 2.0
                ax.annotate("=", xy=(mid, ymax + pad), ha="center", va="bottom",
                           fontsize=16, fontweight="bold", color="#606060",
                           annotation_clip=False)
        if any(is_noop):
            any_noop_rendered = True
            ax.set_ylim(top=ymax + pad * 2.5)

        ax.set_xticks(positions)
        ax.set_xticklabels(present, rotation=30, ha="right")
        ax.set_title(f"N={size}")
        ax.set_ylabel("MFLOPS")
        ax.grid(axis="y", alpha=0.3)

    if any_noop_rendered:
        from matplotlib.patches import Patch
        handles = [Patch(facecolor=color, alpha=0.6, label="change expected"),
                  Patch(facecolor=noop_color, alpha=0.6, hatch="//",
                        edgecolor="#606060",
                        label="\"=\" no change expected (unrelated commit)")]
        fig.legend(handles=handles, loc="lower center", ncol=2,
                  bbox_to_anchor=(0.5, -0.02), frameon=False)

    fig.suptitle(title or f"{algo} {prec} {xform}: MFLOPS distribution across versions")
    fig.tight_layout()
    return fig


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+",
                    help="samples CSV file(s) or glob(s); one version per file "
                         "(grouped by the file's `label` provenance key)")
    ap.add_argument("--order", default=None,
                    help="comma-separated version order for the x-axis "
                         "(default: order first seen among the input files)")
    ap.add_argument("--sizes", required=True,
                    help="comma-separated FFT sizes to plot, one panel each")
    ap.add_argument("--algo", default="pffft", help="algo id (default: pffft)")
    ap.add_argument("--prec", default="flt", choices=["flt", "dbl"])
    ap.add_argument("--xform", default="real", choices=["real", "cplx"])
    ap.add_argument("--kind", default="box", choices=["box", "violin"])
    ap.add_argument("--title", default=None)
    ap.add_argument("--noop", default=None,
                    help="comma-separated version labels expected to show NO "
                         "real change from their predecessor on this "
                         "(prec, xform) series (e.g. a double-only commit "
                         "plotted on the flt series); marked with a hatched "
                         "box and a \"=\" over the gap")
    ap.add_argument("--out", default="bench_results/evolution.png",
                    help="output image path (default: %(default)s)")
    args = ap.parse_args()

    versions = load_versions(args.files)
    if not versions:
        sys.exit("no samples data found in the given files")
    order = args.order.split(",") if args.order else None
    sizes = [int(s) for s in args.sizes.split(",")]

    fig = plot_evolution(versions, sizes, algo=args.algo, prec=args.prec,
                        xform=args.xform, order=order, kind=args.kind,
                        title=args.title,
                        noop=set(args.noop.split(",")) if args.noop else None)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
