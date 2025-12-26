#!/usr/bin/env python3
"""
Convert an extracted PDN SPICE netlist (ALSIM-style) into VoltSpot-style FLP + PTrace + Padloc.

Mode: High-Resolution Tiling
- Instead of grouping by Bxx, this script divides the chip into a fine NxN grid (tiles).
- Loads from SPICE are binned into these tiles to preserve the exact spatial distribution.
- This produces a much more accurate IR drop map that matches the Golden distribution.
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


RE_IB_V = re.compile(r"^\s*iB(\d+)_\d+_v\s+(\S+)\s+0\s+([0-9eE+\-\.]+)\s*$")
RE_NODE_XY = re.compile(r"(?:_X_)?n\d+_(\d+)_(\d+)$")
RE_V_SOURCE = re.compile(r"^\s*V\S+\s+((?:_X_)?n(\d+)_\d+_\d+)\s+0\s+([0-9eE+\-\.]+)\s*$", re.IGNORECASE)
RE_PAD_R = re.compile(r"^\s*R\S+\s+n\d+_\d+_\d+\s+(?:_X_)n\d+_\d+_\d+\s+([0-9eE+\-\.]+)\s*$", re.IGNORECASE)


def node_to_xy(node: str) -> Optional[Tuple[int, int]]:
    m = RE_NODE_XY.search(node)
    if not m:
        return None
    return int(m.group(1)), int(m.group(2))


def parse_spice_data(spice_path: str):
    all_loads = []  # List of (x, y, current)
    pads = []       # List of (type, x, y)
    all_x, all_y = [], []
    detected_vdd = 0.0
    detected_padr = 0.0

    with open(spice_path, "r") as f:
        for line in f:
            # Parse loads
            m_ib = RE_IB_V.match(line)
            if m_ib:
                node, i_a = m_ib.group(2), float(m_ib.group(3))
                xy = node_to_xy(node)
                if xy:
                    all_loads.append((xy[0], xy[1], i_a))
                    all_x.append(xy[0]); all_y.append(xy[1])
                continue

            # Parse Pads
            m_v = RE_V_SOURCE.match(line)
            if m_v:
                node_name, val = m_v.group(1), float(m_v.group(3))
                xy = node_to_xy(node_name)
                if xy:
                    pad_type = "V" if val > 0.1 else "G"
                    pads.append((pad_type, xy[0], xy[1]))
                    if pad_type == "V": detected_vdd = max(detected_vdd, val)
                    all_x.append(xy[0]); all_y.append(xy[1])
                continue

            # Detect PadR
            m_r = RE_PAD_R.match(line)
            if m_r: detected_padr = float(m_r.group(1))

    return all_loads, pads, detected_vdd, detected_padr, all_x, all_y


def main() -> int:
    ap = argparse.ArgumentParser(description="Convert SPICE PDN to VoltSpot inputs (Tiling Mode).")
    ap.add_argument("--spice", required=True)
    ap.add_argument("--num_tiles", type=int, default=100, help="Resolution: number of tiles per side (default 100x100)")
    ap.add_argument("--outdir", default="")
    ap.add_argument("--vdd", type=float, default=0.0)
    ap.add_argument("--chip_width_m", type=float, default=10e-3)
    ap.add_argument("--target_grid_rows", type=int, default=73)
    ap.add_argument("--target_grid_cols", type=int, default=73)
    args = ap.parse_args()

    spice_path = Path(args.spice)
    outdir = Path(args.outdir) if args.outdir else spice_path.parent
    
    all_loads, pads, detected_vdd, detected_padr, all_x, all_y = parse_spice_data(str(spice_path))
    if not all_loads: raise SystemExit("No loads found.")

    vdd = args.vdd if args.vdd > 0.0 else (detected_vdd if detected_vdd > 0.0 else 1.8)
    xmin, xmax = min(all_x), max(all_x)
    ymin, ymax = min(all_y), max(all_y)
    dx, dy = max(1, xmax - xmin), max(1, ymax - ymin)
    chip_w = float(args.chip_width_m)
    chip_h = chip_w * (dy / dx)

    # Tiling Logic
    N = args.num_tiles
    tile_w_spice = dx / N
    tile_h_spice = dy / N
    tile_w_m = chip_w / N
    tile_h_m = chip_h / N

    tile_currents = {} # (tx, ty) -> current_sum
    for x, y, i in all_loads:
        tx = min(N-1, int((x - xmin) / tile_w_spice))
        ty = min(N-1, int((y - ymin) / tile_h_spice))
        tile_currents[(tx, ty)] = tile_currents.get((tx, ty), 0.0) + i

    # Generate Output Files
    out_flp = outdir / "gen.flp"
    out_ptrace = outdir / "gen.ptrace"
    out_padloc = outdir / "gen.padloc"
    outdir.mkdir(parents=True, exist_ok=True)

    # Active tiles only to keep file size sane
    active_tiles = sorted(tile_currents.keys())
    unit_names = [f"T_{tx}_{ty}" for tx, ty in active_tiles]

    # Write FLP
    with open(out_flp, "w") as f:
        f.write("# <unit-name>\t<width>\t<height>\t<left-x>\t<bottom-y>\n")
        for (tx, ty), name in zip(active_tiles, unit_names):
            lx = tx * tile_w_m
            by = ty * tile_h_m
            f.write(f"{name}\t{tile_w_m:.12g}\t{tile_h_m:.12g}\t{lx:.12g}\t{by:.12g}\n")

    # Write PTrace
    with open(out_ptrace, "w") as f:
        f.write("\t".join(unit_names) + "\t\n")
        f.write("\t".join([f"{vdd * tile_currents[t]:.12g}" for t in active_tiles]) + "\t\n")

    # Write Padloc
    mapx_grid = lambda x: int(round((x - xmin) / dx * (args.target_grid_cols - 1)))
    mapy_grid = lambda y: int(round((y - ymin) / dy * (args.target_grid_rows - 1)))
    with open(out_padloc, "w") as f:
        f.write("# <V/G>\t<x_idx>\t<y_idx>\n")
        for ptype, px, py in pads:
            f.write(f"{ptype}\t{mapx_grid(px)}\t{mapy_grid(py)}\n")

    print(f"Tiling completed: {N}x{N} grid, {len(active_tiles)} active units.")
    print(f"Wrote: {out_flp.name}, {out_ptrace.name}, {out_padloc.name} to {outdir}")
    print(f"\nRECOMMENDED pdn.config:\n  -vdd {vdd:g}\n  -PDN_padR {detected_padr if detected_padr > 0 else 0.25:g}\n  -PDN_padconfig 0")
    return 0

if __name__ == "__main__":
    main()
