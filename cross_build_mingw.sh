#!/bin/bash
# Cross-compiles pffft for Windows on Debian or Ubuntu, which need the
# packages gcc-mingw-w64 and g++-mingw-w64. Builds in build_<target>_<name>
# in the current directory.

usage="usage: $0 <w32|w64> <name> <other cmake options>"
case $1 in
  w32) toolchain=mingw-w32-i686.cmake ;;
  w64) toolchain=mingw-w64-x64_64.cmake ;;
  *) echo "$usage"; exit 1 ;;
esac
[ -n "$2" ] || { echo "$usage"; exit 1; }
set -e
src=$(cd "$(dirname "$0")" && pwd)
build=build_$1_$2
install=pffft_bin-$1_$2
shift 2

rm -rf "$build"
cmake -S "$src" -B "$build" -DCMAKE_TOOLCHAIN_FILE="$src/$toolchain" -DCMAKE_INSTALL_PREFIX="$install" "$@"
cmake --build "$build"
