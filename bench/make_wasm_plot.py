#!/usr/bin/env python3
"""Render the WASM-backend evolution plot: last pre-WASM commit (2971ffa,
NEON-emulation Emscripten path) vs. the dedicated WASM SIMD backend commit
(e836207), run under Node via emcc. No noop annotations: this pair isolates
exactly one commit, and it is expected to change WASM throughput by design
(that's the point of the commit) -- there is no "unrelated precision" to
cross-check against here.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "bench"))
from plot_evolution import load_versions, plot_evolution  # noqa: E402

SERIES_DIR = REPO / ".perf" / "series"
OUT_DIR = REPO / "bench_results"
SIZES = [256, 1024, 4096]
ORDER = ["pre-wasm-neon-emu", "wasm-dedicated"]


def render(prec: str, xform: str) -> None:
    pattern = str(SERIES_DIR / f"*-wasm-node-{prec}.csv")
    versions = load_versions([pattern])
    missing = [lab for lab in ORDER if lab not in versions]
    if missing:
        print(f"  [wasm/{prec}/{xform}] skipping, missing labels: {missing}")
        return
    prec_name = {"flt": "float", "dbl": "double"}[prec]
    xform_name = {"real": "real", "cplx": "complex"}[xform]
    fig = plot_evolution(
        versions, SIZES, algo="pffft", prec=prec, xform=xform,
        order=ORDER, kind="box",
        title=(f"pffft {prec_name} {xform_name} MFLOPS: NEON-emulation vs. "
              f"dedicated WASM SIMD backend (Node + emcc)"))
    out = OUT_DIR / f"evolution-wasm-{prec}-{xform}.png"
    fig.savefig(out, dpi=150)
    print(f"  wrote {out}")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for prec in ("flt", "dbl"):
        for xform in ("real", "cplx"):
            render(prec, xform)


if __name__ == "__main__":
    main()
