#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage: ./run.sh testcaseX [vdd]

Example:
  ./run.sh testcase1

This runs:
  - pdn_load_example on benchmarks/<testcase>/gen.flp/ptrace/padloc/config
  - visualize_pdn_loads.py
  - visualize_ir_drop.py
EOF
}

if [[ $# -lt 1 || $# -gt 2 ]]; then
  usage
  exit 1
fi

tc="$1"
vdd="${2:-1.8}"
bench_dir="benchmarks/${tc}"

flp="${bench_dir}/gen.flp"
cfg="${bench_dir}/pdn.config"
ptrace="${bench_dir}/gen.ptrace"
padloc="${bench_dir}/gen.padloc"
outdir="${bench_dir}/out"
mlcf="${bench_dir}/gen.mlcf"
branches_csv="${outdir}/pdn_with_loads_branches.csv"
nodes_csv="${outdir}/pdn_with_loads_branches_nodes.csv"
gridir="${outdir}/ir_drop.gridIR"

if [[ ! -f "$flp" || ! -f "$cfg" ]]; then
  echo "Missing flp or config under ${bench_dir}" >&2
  exit 1
fi

mkdir -p "$outdir"

echo "=== Running pdn_load_example for ${tc} ==="
./build/pdn_load_example \
  --flp "$flp" \
  --config "$cfg" \
  --ptrace "$ptrace" \
  --padloc "$padloc" \
  --outdir "$outdir" \
  --mlcf "$mlcf" \
  --vdd "$vdd" --gnd 0

echo "=== Visualizing PDN loads ==="
python3 tools/visualize_pdn_loads.py \
  --branches "$branches_csv" \
  --nodes "$nodes_csv" \
  --out "${outdir}/pdn_visualization.png"

echo "=== Visualizing IR drop ==="
python3 tools/visualize_ir_drop.py \
  --gridir "$gridir" \
  --out "${outdir}/ir_drop_visualization.png"

echo "Done. Outputs in ${outdir}"

