#!/usr/bin/env python3
"""
Convert an extracted PDN SPICE netlist (ALSIM-style) into VoltSpot-style FLP + PTrace + Padloc.

Target use-case (goal #1):
- You have a benchmark netlist like `benchmarks/testcase1/ibmpg1.spice`
- Loads are modeled as paired current sources iBxx_*_v / iBxx_*_g.
- Pads are modeled as voltage sources V... connecting to _X_ nodes.
- Pad resistance is detected from R... connecting regular nodes to _X_ nodes.

This script generates:
- gen.flp: one rectangular unit per Bxx (bounding box of its load nodes).
- gen.ptrace: P(Bxx) = VDD * sum(I(Bxx)).
- gen.padloc: V/G coordinates mapped to virtual grid indices.
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


@dataclass
class BStats:
    i_sum_a: float
    xmin: int
    ymin: int
    xmax: int
    ymax: int

    @classmethod
    def from_xy_i(cls, x: int, y: int, i: float) -> "BStats":
        return cls(i_sum_a=i, xmin=x, ymin=y, xmax=x, ymax=y)

    def add(self, x: int, y: int, i: float) -> None:
        self.i_sum_a += i
        self.xmin = min(self.xmin, x)
        self.ymin = min(self.ymin, y)
        self.xmax = max(self.xmax, x)
        self.ymax = max(self.ymax, y)


def parse_spice_data(
    spice_path: str,
) -> Tuple[Dict[int, BStats], Dict[int, List[Tuple[int, int]]], List[Tuple[str, int, int]], float, float, List[int], List[int]]:
    stats: Dict[int, BStats] = {}
    points_by_b: Dict[int, List[Tuple[int, int]]] = {}
    pads: List[Tuple[str, int, int]] = []
    all_x: List[int] = []
    all_y: List[int] = []
    detected_vdd = 0.0
    detected_padr = 0.0

    with open(spice_path, "r") as f:
        for line in f:
            # Parse loads
            m_ib = RE_IB_V.match(line)
            if m_ib:
                b = int(m_ib.group(1))
                node = m_ib.group(2)
                i_a = float(m_ib.group(3))
                xy = node_to_xy(node)
                if xy:
                    x, y = xy
                    if b not in stats:
                        stats[b] = BStats.from_xy_i(x, y, i_a)
                    else:
                        stats[b].add(x, y, i_a)
                    points_by_b.setdefault(b, []).append((x, y))
                    all_x.append(x)
                    all_y.append(y)
                continue

            # Parse Pads (Voltage Sources)
            m_v = RE_V_SOURCE.match(line)
            if m_v:
                node_name = m_v.group(1)
                net_id = int(m_v.group(2))
                val = float(m_v.group(3))
                xy = node_to_xy(node_name)
                if xy:
                    x, y = xy
                    # val > 0.1 is VDD, else GND
                    pad_type = "V" if val > 0.1 else "G"
                    pads.append((pad_type, x, y))
                    if pad_type == "V":
                        detected_vdd = max(detected_vdd, val)
                    all_x.append(x)
                    all_y.append(y)
                continue

            # Detect PadR
            m_r = RE_PAD_R.match(line)
            if m_r:
                detected_padr = float(m_r.group(1))

    return stats, points_by_b, pads, detected_vdd, detected_padr, all_x, all_y


def draw_blocks_png(
    *,
    out_png: Path,
    stats: Dict[int, BStats],
    points_by_b: Dict[int, List[Tuple[int, int]]],
    mapx,
    mapy,
    chip_w: float,
    chip_h: float,
    stride: int,
) -> None:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 8 * (chip_h / chip_w)))
    ax.set_title("Bxx block partition + load points")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_xlim(0.0, chip_w)
    ax.set_ylim(0.0, chip_h)
    ax.set_aspect("equal", adjustable="box")

    units = sorted(stats.keys())
    for b in units:
        bb = stats[b]
        lx, by = mapx(bb.xmin), mapy(bb.ymin)
        rx, ty = mapx(bb.xmax), mapy(bb.ymax)
        w, h = max(1e-12, rx - lx), max(1e-12, ty - by)
        ax.add_patch(mpatches.Rectangle((lx, by), w, h, fill=False, linewidth=1.0, alpha=0.8))
        ax.text(lx + 0.5 * w, by + 0.5 * h, f"B{b}", ha="center", va="center", fontsize=8)
        pts = points_by_b.get(b, [])[::max(1, stride)]
        if pts:
            ax.scatter([mapx(x) for x, _ in pts], [mapy(y) for _, y in pts], s=2, alpha=0.1, rasterized=True)

    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description="Convert SPICE PDN to VoltSpot inputs.")
    ap.add_argument("--spice", required=True)
    ap.add_argument("--outdir", default="")
    ap.add_argument("--out_flp", default="")
    ap.add_argument("--out_ptrace", default="")
    ap.add_argument("--out_padloc", default="")
    ap.add_argument("--out_blocks_png", default="")
    ap.add_argument("--blocks_stride", type=int, default=10)
    ap.add_argument("--target_grid_cols", type=int, default=73)
    ap.add_argument("--target_grid_rows", type=int, default=73)
    ap.add_argument("--vdd", type=float, default=0.0)
    ap.add_argument("--chip_width_m", type=float, default=10e-3)
    args = ap.parse_args()

    spice_path = Path(args.spice)
    outdir = Path(args.outdir) if args.outdir else spice_path.parent
    out_flp = Path(args.out_flp) if args.out_flp else (outdir / "gen.flp")
    out_ptrace = Path(args.out_ptrace) if args.out_ptrace else (outdir / "gen.ptrace")
    out_padloc = Path(args.out_padloc) if args.out_padloc else (outdir / "gen.padloc")
    out_png = Path(args.out_blocks_png) if args.out_blocks_png else None

    stats, points_by_b, pads, detected_vdd, detected_padr, all_x, all_y = parse_spice_data(str(spice_path))
    if not stats:
        raise SystemExit("No loads found in SPICE.")

    vdd = args.vdd if args.vdd > 0.0 else (detected_vdd if detected_vdd > 0.0 else 1.8)
    xmin, xmax = min(all_x), max(all_x)
    ymin, ymax = min(all_y), max(all_y)
    dx, dy = max(1, xmax - xmin), max(1, ymax - ymin)
    chip_w = float(args.chip_width_m)
    chip_h = chip_w * (dy / dx)

    mapx = lambda x: (x - xmin) / dx * chip_w
    mapy = lambda y: (y - ymin) / dy * chip_h
    mapx_grid = lambda x: int(round((x - xmin) / dx * (args.target_grid_cols - 1)))
    mapy_grid = lambda y: int(round((y - ymin) / dy * (args.target_grid_rows - 1)))

    # Write FLP
    out_flp.parent.mkdir(parents=True, exist_ok=True)
    with open(out_flp, "w") as f:
        f.write("# <unit-name>\t<width>\t<height>\t<left-x>\t<bottom-y>\n")
        for b in sorted(stats.keys()):
            bb = stats[b]
            lx, by = mapx(bb.xmin), mapy(bb.ymin)
            w, h = max(1e-12, mapx(bb.xmax) - lx), max(1e-12, mapy(bb.ymax) - by)
            f.write(f"B{b}\t{w:.12g}\t{h:.12g}\t{lx:.12g}\t{by:.12g}\n")

    # Write PTrace
    with open(out_ptrace, "w") as f:
        units = sorted(stats.keys())
        f.write("\t".join([f"B{b}" for b in units]) + "\t\n")
        f.write("\t".join([f"{vdd * stats[b].i_sum_a:.12g}" for b in units]) + "\t\n")

    # Write Padloc
    with open(out_padloc, "w") as f:
        f.write("# <V/G>\t<x_idx>\t<y_idx>\n")
        for ptype, px, py in pads:
            f.write(f"{ptype}\t{mapx_grid(px)}\t{mapy_grid(py)}\n")

    if out_png:
        draw_blocks_png(out_png=out_png, stats=stats, points_by_b=points_by_b, mapx=mapx, mapy=mapy, chip_w=chip_w, chip_h=chip_h, stride=args.blocks_stride)

    print(f"Wrote FLP: {out_flp}\nWrote PTrace: {out_ptrace}\nWrote Padloc: {out_padloc}")
    print(f"\n" + "=" * 40 + "\nRECOMMENDED pdn.config CHANGES:\n")
    print(f"  -vdd               {vdd:g}\n  -PDN_padR          {detected_padr if detected_padr > 0 else 0.25:g}")
    print(f"  -PDN_padconfig     0\n  -padloc_file_in    {out_padloc.name}\n" + "=" * 40)
    return 0

if __name__ == "__main__":
    main()
