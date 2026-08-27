#!/usr/bin/env python3
"""Render the perf-branch evolution plots (master -> tip) for docs.

Reads .perf/series/<label>-<platform>-<prec>.csv (written by
collect_series_overlay.py / collect_series_wasm.py) and renders one
box-and-whiskers evolution figure per (platform, prec, xform), with
"=" annotations over any step where the underlying commit's own diff
gives us NO reason to expect a change for that (platform, prec, xform,
size) combination -- e.g. a NEON-double-only commit plotted on the flt
series, or the N=32-only radix-8 codelet plotted at N=256/1024/4096.

See ANALYSIS below for exactly which commit touches what; each entry is
grounded in `git show --stat`/full diff of the actual commit, not guessed.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "bench"))
from plot_evolution import load_versions, plot_evolution  # noqa: E402

SERIES_DIR = REPO / ".perf" / "series"
OUT_DIR = REPO / "bench_results"
SIZES = [32, 256, 1024, 4096]

ORDER = [
    "00-master",
    "01-fma-f32-neon",
    "02-vld1q-neon",
    "03-vtranspose4-armv7-aarch64",
    "04-neon-double-rewrite",
    "05-vmsub-all-backends",
    "06-radix3-ido3-real",
    "07-radix4-ido3-real",
    "08-interleave2-neon-double",
    "09-vtranspose4-aarch64-trn-double",
    "10-radix8-codelet-n32-cplx",
    "11-x86-fma-avx-dbl-sse-flt",
    "12-readme-wasm-tip",
]

# platform -> prec -> xform -> set(labels) expected to show NO real change
# vs. their predecessor, for ALL sizes in SIZES unless overridden below.
# Grounded in each commit's own diff (see bench/PERF_TESTING.md workflow;
# commit shas in bench/collect_series_overlay.py's caller,
# bench/run_evolution_native.sh):
#
#  01 pf_neon_float.h only (vfmaq_f32 FMA)            -> ARM flt only
#  02 pf_neon_float.h only (vld1q_f32 loads)          -> ARM flt only
#  03 pf_neon_float.h only (VTRANSPOSE4: ARMv7 asm +
#     AArch64 trn1/trn2, both new branches)           -> ARM flt only
#  04 pf_neon_double.h + …_from_avx.h (double rewrite) -> ARM dbl only
#  05 VMSUB added to every backend, VCPLXMUL/CONJ
#     rewritten in pf_double.h/pf_float.h (generic)   -> real change everywhere
#  06 pffft_priv_impl.h radix-3 real ido==3 fast path -> real xform only
#  07 pffft_priv_impl.h radf4/radb4 real ido==3       -> real xform only
#  08 pf_neon_double.h INTERLEAVE2/UNINTERLEAVE2       -> ARM dbl only
#  09 pf_neon_double.h VTRANSPOSE4 trn1/trn2; the 1-line
#     pf_neon_float.h hunk is a comment addition only  -> ARM dbl only
#  10 pffft_priv_impl.h radix-8 codelet, gated on
#     Ncvec==8 i.e. N=32 for a complex transform       -> cplx xform, N=32 only
#  11 pf_avx_double.h + pf_sse1_float.h (x86 FMA)      -> x86 (local) only
#  12 README.md only, zero src/ diff                   -> no change anywhere
ALL = {"real", "cplx"}


def uniform(local_flt, local_dbl, ios_flt, ios_dbl):
    """Build a {platform:{prec:{xform:noop_bool}}} entry, same for both xforms."""
    def cell(noop):
        return {"real": noop, "cplx": noop}
    return {
        "local": {"flt": cell(local_flt), "dbl": cell(local_dbl)},
        "ios": {"flt": cell(ios_flt), "dbl": cell(ios_dbl)},
    }


NOOP = {
    "01-fma-f32-neon": uniform(True, True, False, True),
    "02-vld1q-neon": uniform(True, True, False, True),
    "03-vtranspose4-armv7-aarch64": uniform(True, True, False, True),
    "04-neon-double-rewrite": uniform(True, True, True, False),
    "05-vmsub-all-backends": uniform(False, False, False, False),
    "06-radix3-ido3-real": {
        "local": {"flt": {"real": False, "cplx": True}, "dbl": {"real": False, "cplx": True}},
        "ios": {"flt": {"real": False, "cplx": True}, "dbl": {"real": False, "cplx": True}},
    },
    "07-radix4-ido3-real": {
        "local": {"flt": {"real": False, "cplx": True}, "dbl": {"real": False, "cplx": True}},
        "ios": {"flt": {"real": False, "cplx": True}, "dbl": {"real": False, "cplx": True}},
    },
    "08-interleave2-neon-double": uniform(True, True, True, False),
    "09-vtranspose4-aarch64-trn-double": uniform(True, True, True, False),
    # cplx-only, and only at N=32 -- per-size override applied below.
    "10-radix8-codelet-n32-cplx": uniform(True, True, True, True),
    "11-x86-fma-avx-dbl-sse-flt": uniform(False, False, True, True),
    "12-readme-wasm-tip": uniform(True, True, True, True),
}


def noop_map_for(platform: str, prec: str, xform: str) -> dict[int, set[str]]:
    """{size: {noop labels}} for one (platform, prec, xform) evolution plot."""
    labels = {lab for lab, table in NOOP.items()
             if table[platform][prec][xform]}
    m = {size: set(labels) for size in SIZES}
    if xform == "cplx":
        # 10's noop status only holds away from N=32, where Ncvec==8 fires.
        m[32] = m[32] - {"10-radix8-codelet-n32-cplx"}
    return m


def render(platform: str, prec: str, xform: str) -> None:
    pattern = str(SERIES_DIR / f"*-{platform}-{prec}.csv")
    versions = load_versions([pattern])
    present = [lab for lab in ORDER if lab in versions]
    missing = [lab for lab in ORDER if lab not in versions]
    if missing:
        print(f"  [{platform}/{prec}/{xform}] skipping, missing labels: {missing}")
        return
    noop = noop_map_for(platform, prec, xform)
    platform_name = {"local": "Mac (M2, local)", "ios": "iPhone 11 Pro (A13)"}[platform]
    prec_name = {"flt": "float", "dbl": "double"}[prec]
    xform_name = {"real": "real", "cplx": "complex"}[xform]
    fig = plot_evolution(
        versions, SIZES, algo="pffft", prec=prec, xform=xform,
        order=present, kind="box", noop=noop,
        title=(f"pffft {prec_name} {xform_name} MFLOPS: master -> perf tip "
              f"({platform_name})"))
    out = OUT_DIR / f"evolution-{platform}-{prec}-{xform}.png"
    fig.savefig(out, dpi=150)
    print(f"  wrote {out}")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for platform in ("local", "ios"):
        for prec in ("flt", "dbl"):
            for xform in ("real", "cplx"):
                render(platform, prec, xform)


if __name__ == "__main__":
    main()
