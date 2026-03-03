#!/bin/bash
# Generate info.json for benchmark results directory
set -e

if [ -z "$1" ]; then
  echo "Usage: $0 <output-dir>" >&2
  exit 1
fi

outdir="$1"
mkdir -p "$outdir"

# Detect CPU
if [ "$(uname)" = "Darwin" ]; then
  cpu=$(sysctl -n machdep.cpu.brand_string 2>/dev/null || echo "unknown")
else
  cpu=$(grep -m1 'model name' /proc/cpuinfo 2>/dev/null | cut -d: -f2 | xargs || echo "unknown")
fi

arch=$(uname -m)
os="$(uname -s) $(uname -r)"

# Detect compiler
cc="${CC:-cc}"
compiler=$($cc --version 2>&1 | head -1 || echo "unknown")

date=$(date +%Y-%m-%d)

cat > "$outdir/info.json" <<EOF
{
  "cpu": "$cpu",
  "arch": "$arch",
  "os": "$os",
  "compiler": "$compiler",
  "date": "$date"
}
EOF

echo "Wrote $outdir/info.json"
