#!/bin/bash
# Checks whether a riscv64 build of pffft has vector code. The RVV backend is
# written with GCC vector types, which a compiler may lower to scalar code;
# vsetvli/vsetivli start every vector operation.
#
# usage: check-rvv-code.sh <libpffft.a> <some|none>
# environment:
#   OBJDUMP  an objdump for riscv64, default riscv64-linux-gnu-objdump
set -euo pipefail
lib=$1
expect=$2
objdump=${OBJDUMP:-riscv64-linux-gnu-objdump}

disasm=$("$objdump" -d "$lib") \
  || { echo "::error::$objdump could not disassemble $lib"; exit 1; }
grep -q '<pffft_new_setup>:' <<< "$disasm" \
  || { echo "::error::no pffft code in the disassembly of $lib"; exit 1; }
n=$(grep -cE '\svseti?vli\s' <<< "$disasm" || true)
echo "vsetvli/vsetivli in $lib: $n"
case $expect in
  some) [ "$n" -gt 0 ] || { echo "::error::expected vector instructions in $lib, found none"; exit 1; } ;;
  none) [ "$n" -eq 0 ] || { echo "::error::expected no vector instructions in $lib, found $n"; exit 1; } ;;
  *) echo "usage: $0 <libpffft.a> <some|none>" >&2; exit 2 ;;
esac
