#!/usr/bin/env python3
"""Plain-language before/after MFLOPS summary, for a README, not an audit log.

Renders exactly two bars per panel -- baseline vs final -- with IQR whiskers
and a shaded "noise floor" band around the baseline bar. A reader doesn't
need to know what Mann-Whitney U is: if the final bar's whisker clears the
shaded band, the change is bigger than measurement noise; if it doesn't,
it isn't (yet) distinguishable from noise, whatever the raw percentage says.

The noise floor itself should come from a genuine same-code control
comparison (see PERF_TESTING.md): two measurements of the SAME binary,
at (prec, xform, size) combinations no candidate change could plausibly
touch, run far enough apart in the session to capture real-world drift.
Pass it explicitly with --noise-floor-pct once you've measured it for
your platform; there is no universally-correct default.

Usage
-----
    ./bench/plot_summary.py baseline.csv final.csv \\
        --baseline-label "master" --final-label "bench-update" \\
        --platform adb --sizes 256,1024,4096 --prec flt --xform real \\
        --noise-floor-pct 0.5 \\
        --out bench_results/summary-adb-flt-real.png
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from samples import read_samples  # noqa: E402
from stats import compare, trim   # noqa: E402


def select(rows, algo, prec, xform, size):
    return [s.mflops for s in rows
            if s.algo == algo and s.prec == prec and s.xform == xform and s.size == size]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("baseline_file")
    ap.add_argument("final_file")
    ap.add_argument("--baseline-label", default="baseline")
    ap.add_argument("--final-label", default="final")
    ap.add_argument("--platform", default="")
    ap.add_argument("--sizes", required=True, help="comma-separated FFT sizes, one panel each")
    ap.add_argument("--algo", default="pffft")
    ap.add_argument("--prec", default="flt", choices=["flt", "dbl"])
    ap.add_argument("--xform", default="real", choices=["real", "cplx"])
    ap.add_argument("--noise-floor-pct", type=float, required=True,
                    help="measurement noise floor as a percent of baseline "
                    "MFLOPS (from a genuine same-code control comparison; "
                    "see PERF_TESTING.md); shaded band drawn at +/- this")
    ap.add_argument("--out", default="bench_results/summary.png")
    ap.add_argument("--caption", default=None,
                    help="small caveat/footnote text rendered below the panels "
                    "(e.g. noting a baseline substitution)")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _, base_rows = read_samples(args.baseline_file)
    _, final_rows = read_samples(args.final_file)
    sizes = [int(s) for s in args.sizes.split(",")]

    fig, axes = plt.subplots(1, len(sizes), figsize=(3.2 * len(sizes) + 1, 4.2))
    if len(sizes) == 1:
        axes = [axes]

    for ax, size in zip(axes, sizes):
        bvals = np.array(select(base_rows, args.algo, args.prec, args.xform, size))
        fvals = np.array(select(final_rows, args.algo, args.prec, args.xform, size))
        if len(bvals) == 0 or len(fvals) == 0:
            ax.set_title(f"N={size}\n(no data)")
            continue
        bt, ft = trim(bvals), trim(fvals)
        bmed, fmed = np.median(bt), np.median(ft)
        biqr = (np.percentile(bt, 75) - np.percentile(bt, 25)) / 2
        fiqr = (np.percentile(ft, 75) - np.percentile(ft, 25)) / 2

        band = bmed * args.noise_floor_pct / 100.0
        ax.axhspan(bmed - band, bmed + band, color="0.85", zorder=0,
                  label=f"noise floor (+/-{args.noise_floor_pct:g}%)")

        c = compare(bvals.tolist(), fvals.tolist())
        pct = (fmed - bmed) / bmed * 100
        clears = abs(pct) > args.noise_floor_pct
        color = "#2ca02c" if (clears and pct > 0) else ("#d62728" if (clears and pct < 0) else "#7f7f7f")

        ax.bar([0, 1], [bmed, fmed], yerr=[biqr, fiqr], capsize=6,
              color=["#9ecae1", color], edgecolor="black", width=0.6, zorder=2)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([args.baseline_label, args.final_label], rotation=0, ha="center")
        verdict = "faster" if (clears and pct > 0) else ("slower" if (clears and pct < 0) else "within noise")
        ax.set_title(f"N={size}\n{pct:+.1f}% ({verdict})")
        if ax is axes[0]:
            ax.set_ylabel("MFLOPS")

    title = f"{args.platform + ': ' if args.platform else ''}{args.algo} {args.prec} {args.xform} -- before vs after"
    fig.suptitle(title)
    axes[0].legend(loc="lower left", fontsize=8, framealpha=0.9)
    if args.caption:
        fig.tight_layout(rect=(0, 0.08, 1, 1))
        fig.text(0.5, 0.01, args.caption, ha="center", va="bottom",
                 fontsize=8, style="italic", color="0.3", wrap=True)
    else:
        fig.tight_layout()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
