# PDN SPICE tools

## `spice_netviz.py`

Visualize per-layer resistor networks from an extracted SPICE PDN netlist (e.g. `PDN_research/PDN_Sensitivity_Analysis/spice/ibmpg1.spice`) and export edges to CSV.

### List layers

```bash
python3 PDN_research/PDN_Sensitivity_Analysis/tools/spice_netviz.py --netlist PDN_research/PDN_Sensitivity_Analysis/spice/ibmpg1.spice --list-layers
```

### Render one layer (SVG)

```bash
python3 PDN_research/PDN_Sensitivity_Analysis/tools/spice_netviz.py \
  --netlist PDN_research/PDN_Sensitivity_Analysis/spice/ibmpg1.spice \
  --layer "M5,VDD" \
  --out PDN_research/PDN_Sensitivity_Analysis/out_svg/m5_vdd.svg
```

### Crop to a region (faster)

`--bbox xmin xmax ymin ymax` keeps only edges whose endpoints both lie inside the window.

```bash
python3 PDN_research/PDN_Sensitivity_Analysis/tools/spice_netviz.py \
  --netlist PDN_research/PDN_Sensitivity_Analysis/spice/ibmpg1.spice \
  --layer "M6,GND" \
  --bbox 0 10000 0 10000 \
  --out PDN_research/PDN_Sensitivity_Analysis/out_svg/m6_gnd_crop.svg
```

### Render all layers + export CSV

```bash
python3 PDN_research/PDN_Sensitivity_Analysis/tools/spice_netviz.py \
  --netlist PDN_research/PDN_Sensitivity_Analysis/spice/ibmpg1.spice \
  --all \
  --outdir PDN_research/PDN_Sensitivity_Analysis/out_svg/all \
  --export-csv PDN_research/PDN_Sensitivity_Analysis/out_svg/csv
```

### 3D visualization (Matplotlib)

Exports a simple 3D line plot (layers separated along Z; vias as vertical lines) as a PNG.

```bash
python3 PDN_research/PDN_Sensitivity_Analysis/tools/spice_netviz.py \
  --netlist PDN_research/PDN_Sensitivity_Analysis/spice/ibmpg1.spice \
  --out-3d-png PDN_research/PDN_Sensitivity_Analysis/out_3d/pdn.png
```

### 3D interactive view (Matplotlib)

Opens an interactive window you can rotate with the mouse:

```bash
python3 PDN_research/PDN_Sensitivity_Analysis/tools/spice_netviz.py \
  --netlist PDN_research/PDN_Sensitivity_Analysis/spice/ibmpg1.spice \
  --show-3d \
  --3d-stride 5
```

Optional:

- `--z-step 2000` adjusts layer separation (default `2000`).
- `--no-vias` omits via connections.
- `--3d-stride 5` plots every 5th metal segment (much faster for big nets).
- `--vias-stride 25 --vias-alpha 0.08` reduces the number/visibility of via lines (less clutter).

## Pixel model export (VDD)

Exports a simple pixelized VDD model from the extracted SPICE netlist:

- `rx_m5.csv`, `ry_m5.csv`: per-pixel Rx/Ry for `M5,VDD`
- `rx_m6.csv`, `ry_m6.csv`: per-pixel Rx/Ry for `M6,VDD`
- `via_count.csv`, `rz_via.csv`: per-pixel via count and `Rz = rz_per_via / via_count` between M5 and M6
- `meta.json`: pitch/origin/grid size and basic stats

```bash
python3 PDN_research/PDN_Sensitivity_Analysis/tools/spice_netviz.py \
  --netlist PDN_research/PDN_Sensitivity_Analysis/spice/ibmpg1.spice \
  --export-pixel-vdd PDN_research/PDN_Sensitivity_Analysis/out_pixel/vdd \
  --target-pixels 64 \
  --rz-per-via 1e-4
```

### Pixel model export (VDD, merged single-layer)

Merges `M5,VDD` and `M6,VDD` into a single in-plane sheet (parallel conductance sum) and exports:

- `rx.csv`, `ry.csv`, `rz.csv`, `via_count.csv`, `meta.json`

```bash
python3 PDN_research/PDN_Sensitivity_Analysis/tools/spice_netviz.py \
  --netlist PDN_research/PDN_Sensitivity_Analysis/spice/ibmpg1.spice \
  --export-pixel-vdd-merged PDN_research/PDN_Sensitivity_Analysis/out_pixel/vdd_merged \
  --target-pixels 64 \
  --rz-per-via 1e-4
```

### Pixel heatmaps

```bash
python3 PDN_research/PDN_Sensitivity_Analysis/tools/spice_netviz.py \
  --plot-pixel PDN_research/PDN_Sensitivity_Analysis/out_pixel/vdd_merged
```

## Visualize virtual grid (branches.csv)

If you built a virtual grid via the C++ examples and got a `branches.csv` (e.g. `PDN_research/PDN_Sensitivity_Analysis/out_voltspot/branches.csv`),
you can visualize it with Matplotlib:

```bash
python3 PDN_research/PDN_Sensitivity_Analysis/tools/spice_netviz.py \
  --plot-branches PDN_research/PDN_Sensitivity_Analysis/out_voltspot/branches.csv \
  --show-branches \
  --branches-stride 1 \
  --branches-alpha 0.25
```

### Notes

- This assumes node names follow the common extracted-PDN convention `n<net>_<x>_<y>` (and `_X_n<net>_<x>_<y>`).
- SVG edge colors map resistor values (log scale) from low (blue) to high (red).
- For very large nets, use `--bbox` and/or `--max-edges` to keep output sizes reasonable.

## `spice_to_flp_ptrace.py`

Generate VoltSpot-style `*.flp` + `*.ptrace` from an extracted PDN SPICE netlist by grouping loads into `Bxx` units.

This expects current sources named like:

- `iB33_0_v  <VDD_node> 0 <I>`
- `iB33_0_g  0 <GND_node> <I>` (ignored for power summation to avoid double counting)

Power is computed as: `P(Bxx) = VDD * sum(I(Bxx))`.

Example (testcase1):

```bash
python3 tools/spice_to_flp_ptrace.py \
  --spice benchmarks/testcase1/ibmpg1.spice \
  --vdd 1.8 \
  --out_blocks_png benchmarks/testcase1/blocks.png
```

By default, outputs are written next to the input netlist:

- `benchmarks/testcase1/gen.flp`
- `benchmarks/testcase1/gen.ptrace`
