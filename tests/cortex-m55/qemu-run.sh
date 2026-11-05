#!/bin/sh
# usage: qemu-run.sh <program.elf> [args...]
# Runs a semihosting program on QEMU's mps3-an547 (Cortex-M55); the exit
# code is the program's.
elf=$1; shift
cfg="enable=on,target=native,arg=$(basename "$elf")"
for a in "$@"; do cfg="$cfg,arg=$(printf '%s' "$a" | sed 's/,/,,/g')"; done
exec qemu-system-arm -M mps3-an547 -nographic -monitor none -serial none \
  -semihosting-config "$cfg" -kernel "$elf"
