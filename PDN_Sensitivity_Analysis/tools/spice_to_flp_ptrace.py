#!/usr/bin/env python3
"""
Convert an extracted PDN SPICE netlist (ALSIM-style) into VoltSpot-style FLP + PTrace + Padloc.

Mode: Block-Based Extraction (Bxx)
- Groups all iB<num>_... current sources into a single functional unit 'B<num>'.
- Calculates the bounding box of all points within a block to define the FLP unit.
- Strictly follows the VDD net-index filtering logic.
- Visualizes the floorplan in blocks.png.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set


# Regex for strictly following the specification
RE_LAYER = re.compile(r"^\*\s*layer:\s*[^,]+,\s*(\w+)\s*net:\s*(\d+)", re.IGNORECASE)
RE_NODE_XY = re.compile(r"n(\d+)_(\d+)_(\d+)")
RE_IB = re.compile(r"^\s*iB(\d+)\S*\s+(\S+)\s+(\S+)\s+([0-9eE+\-\.]+)", re.IGNORECASE)
RE_V_SOURCE = re.compile(r"^\s*V\S+\s+(\S+)\s+0\s+([0-9eE+\-\.]+)\s*$", re.IGNORECASE)
RE_PAD_R = re.compile(r"^\s*R\S+\s+n\d+_\d+_\d+\s+(?:_X_)n\d+_\d+_\d+\s+([0-9eE+\-\.]+)\s*$", re.IGNORECASE)


def parse_node(node: str) -> Optional[Tuple[int, int, int]]:
    """Returns (net_index, x, y) from node name."""
    m = RE_NODE_XY.search(node)
    if not m:
        return None
    return int(m.group(1)), int(m.group(2)), int(m.group(3))


def parse_spice_data(spice_path: str):
    # block_id -> { 'points': [(x,y)], 'current': total_i }
    blocks_data = {}
    pads = []
    all_x, all_y = [], []
    detected_vdd = 0.0
    detected_padr = 0.0
    
    vdd_nets: Set[int] = set()

    with open(spice_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # 1. Parse Layer/Net definitions
            if line.startswith("*"):
                m_layer = RE_LAYER.match(line)
                if m_layer:
                    net_type = m_layer.group(1).upper()
                    net_idx = int(m_layer.group(2))
                    if "VDD" in net_type:
                        vdd_nets.add(net_idx)
                continue

            # 2. Parse Current Sources (iB...) - Grouped by B<num>
            if line.lower().startswith("ib"):
                m_ib = RE_IB.match(line)
                if m_ib:
                    b_id, node1, node2, val = m_ib.group(1), m_ib.group(2), m_ib.group(3), float(m_ib.group(4))
                    
                    info1 = parse_node(node1)
                    if info1 and info1[0] in vdd_nets and node2 == "0":
                        x, y = info1[1], info1[2]
                        unit_name = f"B{b_id}"
                        if unit_name not in blocks_data:
                            blocks_data[unit_name] = {'pts': [], 'curr': 0.0}
                        blocks_data[unit_name]['pts'].append((x, y))
                        blocks_data[unit_name]['curr'] += val
                        all_x.append(x); all_y.append(y)
                continue

            # 3. Parse Voltage Sources (Pads)
            if line.lower().startswith("v"):
                m_v = RE_V_SOURCE.match(line)
                if m_v:
                    node_name, val = m_v.group(1), float(m_v.group(2))
                    if "_X_" in node_name:
                        info = parse_node(node_name)
                        if info:
                            pad_type = "V" if val > 0.1 else "G"
                            pads.append((pad_type, info[1], info[2]))
                            if pad_type == "V":
                                detected_vdd = max(detected_vdd, val)
                            all_x.append(info[1]); all_y.append(info[2])
                continue

            # 4. Detect Pad Resistance
            if line.lower().startswith("r"):
                if "_X_" in line:
                    m_r = RE_PAD_R.match(line)
                    if m_r:
                        detected_padr = float(m_r.group(1))

    return blocks_data, pads, detected_vdd, detected_padr, all_x, all_y


def draw_blocks_png(blocks_data, pads, xmin, xmax, ymin, ymax, output_path: Path):
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
    except ImportError:
        print("Warning: matplotlib not found, skipping blocks.png generation.")
        return

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_title("Extracted Bxx Blocks & Pads Locations")
    ax.set_xlabel("X (Spice Units)")
    ax.set_ylabel("Y (Spice Units)")

    # Draw Blocks
    for name, data in blocks_data.items():
        pts = data['pts']
        b_xmin = min(p[0] for p in pts)
        b_xmax = max(p[0] for p in pts)
        b_ymin = min(p[1] for p in pts)
        b_ymax = max(p[1] for p in pts)
        
        width = max(1, b_xmax - b_xmin)
        height = max(1, b_ymax - b_ymin)
        
        rect = patches.Rectangle((b_xmin, b_ymin), width, height, 
                                 linewidth=1, edgecolor='blue', facecolor='cyan', alpha=0.3)
        ax.add_patch(rect)
        ax.text(b_xmin + width/2, b_ymin + height/2, name, 
                ha='center', va='center', fontsize=8, color='blue')

    # Draw Pads
    v_pads_x = [p[1] for p in pads if p[0] == "V"]
    v_pads_y = [p[2] for p in pads if p[0] == "V"]
    g_pads_x = [p[1] for p in pads if p[0] == "G"]
    g_pads_y = [p[2] for p in pads if p[0] == "G"]
    
    ax.scatter(v_pads_x, v_pads_y, c='red', s=10, label='VDD Pad', marker='s')
    ax.scatter(g_pads_x, g_pads_y, c='black', s=10, label='GND Pad', marker='x')
    
    ax.legend()
    plt.savefig(output_path)
    plt.close()
    print(f"Generated visualization: {output_path}")


def main() -> int:
    ap = argparse.ArgumentParser(description="Convert SPICE PDN to VoltSpot inputs (Block Mode).")
    ap.add_argument("--spice", required=True)
    ap.add_argument("--outdir", default="")
    ap.add_argument("--vdd", type=float, default=0.0)
    ap.add_argument("--chip_width_m", type=float, default=10e-3)
    ap.add_argument("--target_grid_rows", type=int, default=73)
    ap.add_argument("--target_grid_cols", type=int, default=73)
    args = ap.parse_args()

    spice_path = Path(args.spice)
    outdir = Path(args.outdir) if args.outdir else spice_path.parent
    
    blocks_data, pads, detected_vdd, detected_padr, all_x, all_y = parse_spice_data(str(spice_path))
    if not blocks_data:
        print("Error: No current loads found.")
        return 1

    vdd = args.vdd if args.vdd > 0.0 else (detected_vdd if detected_vdd > 0.0 else 1.8)
    xmin, xmax = min(all_x), max(all_x)
    ymin, ymax = min(all_y), max(all_y)
    dx, dy = max(1, xmax - xmin), max(1, ymax - ymin)
    chip_w = float(args.chip_width_m)
    chip_h = chip_w * (dy / dx)
    
    scale = chip_w / dx

    # Generate Output Files
    out_flp = outdir / "gen.flp"
    out_ptrace = outdir / "gen.ptrace"
    out_padloc = outdir / "gen.padloc"
    out_png = outdir / "blocks.png"
    outdir.mkdir(parents=True, exist_ok=True)

    block_names = sorted(blocks_data.keys(), key=lambda x: int(x[1:]))

    # Write FLP
    with open(out_flp, "w") as f:
        f.write("# <unit-name>\t<width>\t<height>\t<left-x>\t<bottom-y>\n")
        for name in block_names:
            pts = blocks_data[name]['pts']
            b_xmin = min(p[0] for p in pts)
            b_xmax = max(p[0] for p in pts)
            b_ymin = min(p[1] for p in pts)
            b_ymax = max(p[1] for p in pts)
            w_m = max(1e-7, (b_xmax - b_xmin) * scale)
            h_m = max(1e-7, (b_ymax - b_ymin) * scale)
            lx_m = (b_xmin - xmin) * scale
            by_m = (b_ymin - ymin) * scale
            f.write(f"{name}\t{w_m:.12g}\t{h_m:.12g}\t{lx_m:.12g}\t{by_m:.12g}\n")

    # Write PTrace
    with open(out_ptrace, "w") as f:
        f.write("\t".join(block_names) + "\n")
        f.write("\t".join([f"{vdd * blocks_data[name]['curr']:.12g}" for name in block_names]) + "\n")

    # Write Padloc
    mapx_grid = lambda x: int(round((x - xmin) / dx * (args.target_grid_cols - 1)))
    mapy_grid = lambda y: int(round((y - ymin) / dy * (args.target_grid_rows - 1)))
    with open(out_padloc, "w") as f:
        f.write("# <V/G>\t<x_idx>\t<y_idx>\n")
        for ptype, px, py in pads:
            f.write(f"{ptype}\t{mapx_grid(px)}\t{mapy_grid(py)}\n")

    # Visualize
    draw_blocks_png(blocks_data, pads, xmin, xmax, ymin, ymax, out_png)

    print(f"Extraction completed (Block Mode: Bxx).")
    print(f"Total Blocks: {len(block_names)}")
    print(f"Total Current extracted: {sum(d['curr'] for d in blocks_data.values()):.3f} A")
    print(f"Wrote: {out_flp.name}, {out_ptrace.name}, {out_padloc.name} to {outdir}")
    print(f"\nRECOMMENDED pdn.config:\n  -vdd {vdd:g}\n  -PDN_padR {detected_padr if detected_padr > 0 else 0.25:g}\n  -PDN_padconfig 0")
    return 0


if __name__ == "__main__":
    main()
