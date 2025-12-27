#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: ./gen.sh testcaseX [vdd]

Example:
  ./gen.sh testcase1 1.8

Runs spice_to_flp_ptrace.py for the given testcase under benchmarks/.
If vdd is omitted, the script will let the extractor auto-detect or default.
EOF
}

if [[ $# -lt 1 || $# -gt 2 ]]; then
  usage
  exit 1
fi

tc="$1"
vdd_arg=""
if [[ $# -eq 2 ]]; then
  vdd_arg="--vdd $2"
fi

suffix="${tc#testcase}"
if [[ "$suffix" == "$tc" || -z "$suffix" ]]; then
  spice="benchmarks/${tc}/ibmpg1.spice"
else
  spice="benchmarks/${tc}/ibmpg${suffix}.spice"
fi
outdir="benchmarks/${tc}"
blocks_png="${outdir}/blocks.png"

if [[ ! -f "$spice" ]]; then
  echo "SPICE file not found: $spice" >&2
  exit 1
fi

python3 tools/spice_to_flp_ptrace.py \
  --spice "$spice" \
  $vdd_arg \
  --outdir "$outdir"

echo "Done. Generated FLP/PTrace/Padloc under $outdir"

