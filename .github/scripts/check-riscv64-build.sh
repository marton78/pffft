#!/bin/bash
# Checks a riscv64 build of pffft against what a CI job expects.
#
# usage: check-riscv64-build.sh <build dir> <configure log>
# environment:
#   BACKEND  the SIMD backend that cmake must report for float and double;
#            "RVV (scalar code)" for an RVV backend without vector
#            instructions, which cmake must also warn about
#   MARCH    the only -march the build may use; empty for none
#   OBJDUMP  an objdump for riscv64, default riscv64-linux-gnu-objdump
set -euo pipefail
build=$1
log=$2

grep "SIMD backend for PFFFT" "$log"
grep -qF "SIMD backend for PFFFT: float: $BACKEND, double: $BACKEND" "$log" \
  || { echo "::error::expected the $BACKEND backend"; exit 1; }
if [ "$BACKEND" = "RVV (scalar code)" ]; then
  # cmake wraps the lines of a warning
  tr -s ' \n' '  ' < "$log" | grep -q "compiles pffft's RVV backend to scalar code" \
    || { echo "::error::expected a warning about the scalar RVV code"; exit 1; }
fi

db=$build/compile_commands.json
[ -s "$db" ] || { echo "::error::$db is missing or empty"; exit 1; }
# grep exits with 1 when nothing matches, which is fine here
found=$({ grep -o -- '-march=[^ "]*' "$db" || [ $? -eq 1 ]; } | sort -u)
expected=${MARCH:+-march=$MARCH}
echo "-march found: '$found', expected: '$expected'"
[ "$found" = "$expected" ] || { echo "::error::unexpected -march"; exit 1; }

if [ "$BACKEND" = RVV ]; then expect=some; else expect=none; fi
bash "$(dirname "$0")/check-rvv-code.sh" "$build/libpffft.a" $expect
