#!/usr/bin/env bash
# Resume run_evolution_native.sh from step 01 (00-master already collected).
set -euo pipefail
cd "$(dirname "$0")/.."

SIZES="32,256,1024,4096"
RUNS=20
WARMUP=10

STEPS=(
  "01-fma-f32-neon:008b209741f8426a743e465fc47803f57a27d6c7"
  "02-vld1q-neon:d3cf0f9bcbf924fc4429477ab5f15e26c24c869d"
  "03-vtranspose4-armv7-aarch64:e3e6d41344f3b123075efc0630bb4f879386de05"
  "04-neon-double-rewrite:e20a3f0a9bd5f31705dc75c1e3c5143297c1ad04"
  "05-vmsub-all-backends:60a1b8dcd8f3fc0a9dde2be3ec5b3063d8cb9337"
  "06-radix3-ido3-real:cf7d5a090243570e8e077e4294de9b875d9947c5"
  "07-radix4-ido3-real:1e1f5750ea439f0423b1870870555ea875f04245"
  "08-interleave2-neon-double:78825ef778c7caf95822f591e2b365dce4db7e6d"
  "09-vtranspose4-aarch64-trn-double:9f311d1d5c6a85758e44fafae41d61fed8ef8004"
  "10-radix8-codelet-n32-cplx:41683ac1cbf140b8d0e3d3199cd78c4ea5c42179"
  "11-x86-fma-avx-dbl-sse-flt:2971ffae8614ad3599e4312aff35a0f144e37df5"
  "12-readme-wasm-tip:5b220c08f580aed12f96debce3fbb0ff9389e700"
)

for step in "${STEPS[@]}"; do
  label="${step%%:*}"
  ref="${step#*:}"
  echo "=========================================================="
  echo "== $(date '+%H:%M:%S') step $label ($ref)"
  echo "=========================================================="
  python3 bench/collect_series_overlay.py \
    --ref "$ref" --harness-ref bench-update --label "$label" \
    --target local --target ios \
    --sizes "$SIZES" --runs "$RUNS" --warmup-runs "$WARMUP" \
    --prec both --wt .perf/wt-seq --out .perf/series \
    --serial-measure
done

echo "== native evolution sweep complete =="
