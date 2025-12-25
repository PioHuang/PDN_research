#!/usr/bin/env python3
"""
Convert an extracted PDN SPICE netlist (ALSIM-style) into VoltSpot-style FLP + PTrace.

Target use-case (goal #1):
- You have a benchmark netlist like `benchmarks/testcase1/ibmpg1.spice`
- Loads are modeled as paired current sources named like:
    iB33_0_v  <VDD_node> 0  <I>
    iB33_0_g  0 <GND_node> <I>
- Pads exist, but for the C++ flow you can keep PDN_padconfig=1 and skip padloc.

This script generates:
- gen.flp: one rectangular unit per Bxx (bounding box of its *_v current source nodes)
- gen.ptrace: one time step with P(Bxx) = VDD * sum(I(Bxx))  (W)

Notes:
- The geometry in the netlist is usually in integer "database units" in node names, e.g. n3_7130_471.
  This script maps them into meters by preserving aspect ratio and using a chosen chip width in meters.
- This is an approximation to enable VoltSpot-style flows; it will not reproduce the original irregular PDN exactly.
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


RE_IB_V = re.compile(r"^\s*iB(\d+)_\d+_v\s+(\S+)\s+0\s+([0-9eE+\-\.]+)\s*$")
RE_NODE_XY = re.compile(r"(?:_X_)?n\d+_(\d+)_(\d+)$")


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


def parse_spice_currents(
    spice_path: str,
) -> Tuple[Dict[int, BStats], Dict[int, List[Tuple[int, int]]], List[int], List[int]]:
    stats: Dict[int, BStats] = {}
    points_by_b: Dict[int, List[Tuple[int, int]]] = {}
    all_x: List[int] = []
    all_y: List[int] = []

    with open(spice_path, "r") as f:
        for line in f:
            m = RE_IB_V.match(line)
            if not m:
                continue

            b = int(m.group(1))
            node = m.group(2)
            i_a = float(m.group(3))

            xy = node_to_xy(node)
            if xy is None:
                continue
            x, y = xy

            if b not in stats:
                stats[b] = BStats.from_xy_i(x, y, i_a)
            else:
                stats[b].add(x, y, i_a)

            points_by_b.setdefault(b, []).append((x, y))
            all_x.append(x)
            all_y.append(y)

    return stats, points_by_b, all_x, all_y


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
    # Lazy import to avoid requiring matplotlib unless visualization is requested.
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    out_png.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 8 * (chip_h / chip_w)))
    ax.set_title("Bxx block partition (bbox from iBxx_*_v load locations)")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_xlim(0.0, chip_w)
    ax.set_ylim(0.0, chip_h)
    ax.set_aspect("equal", adjustable="box")

    # Draw bboxes + (optionally subsampled) points.
    units = sorted(stats.keys())
    for b in units:
        bb = stats[b]
        lx = mapx(bb.xmin)
        by = mapy(bb.ymin)
        rx = mapx(bb.xmax)
        ty = mapy(bb.ymax)
        w = max(1e-12, rx - lx)
        h = max(1e-12, ty - by)

        rect = mpatches.Rectangle((lx, by), w, h, fill=False, linewidth=1.0, alpha=0.9)
        ax.add_patch(rect)
        ax.text(lx + 0.5 * w, by + 0.5 * h, f"B{b}", ha="center", va="center", fontsize=8, alpha=0.9)

        pts = points_by_b.get(b, [])
        if pts:
            s = max(1, stride)
            pts = pts[::s]
            xs = [mapx(x) for x, _ in pts]
            ys = [mapy(y) for _, y in pts]
            ax.scatter(xs, ys, s=2, alpha=0.15, rasterized=True)

    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Convert ALSIM-style PDN SPICE netlist loads (iBxx_*_v) into VoltSpot FLP + PTrace (Bxx units)."
    )
    ap.add_argument("--spice", required=True, help="Input SPICE netlist (e.g. benchmarks/testcase1/ibmpg1.spice)")
    ap.add_argument(
        "--outdir",
        default="",
        help="Output directory. Default: directory containing --spice (recommended: benchmarks/testcase1).",
    )
    ap.add_argument("--out_flp", default="", help="Output FLP path. Default: <outdir>/gen.flp")
    ap.add_argument("--out_ptrace", default="", help="Output PTrace path. Default: <outdir>/gen.ptrace")
    ap.add_argument("--out_blocks_png", default="", help="Optional: output PNG visualizing Bxx bboxes and load points.")
    ap.add_argument(
        "--blocks_stride",
        type=int,
        default=1,
        help="Point subsampling stride for --out_blocks_png (default: 10, i.e. plot every 10th load point).",
    )
    ap.add_argument("--vdd", type=float, default=1.8, help="Supply voltage used to convert I->P (default: 1.8)")
    ap.add_argument(
        "--chip_width_m",
        type=float,
        default=10e-3,
        help="Chip width in meters used to scale coordinates (default: 10e-3)",
    )
    args = ap.parse_args()

    spice_path = Path(args.spice)
    outdir = Path(args.outdir) if args.outdir else spice_path.parent
    out_flp = Path(args.out_flp) if args.out_flp else (outdir / "gen.flp")
    out_ptrace = Path(args.out_ptrace) if args.out_ptrace else (outdir / "gen.ptrace")
    out_png = Path(args.out_blocks_png) if args.out_blocks_png else None

    stats, points_by_b, all_x, all_y = parse_spice_currents(str(spice_path))
    if not stats:
        raise SystemExit(
            "No matching current sources found. Expected lines like: iB33_0_v <node> 0 <I>\n"
            "If your netlist uses a different naming convention, adjust the regex in this script."
        )

    if not all_x or not all_y:
        raise SystemExit("Parsed current sources but could not extract any (x,y) from node names.")

    xmin, xmax = min(all_x), max(all_x)
    ymin, ymax = min(all_y), max(all_y)
    dx = max(1, xmax - xmin)
    dy = max(1, ymax - ymin)

    chip_w = float(args.chip_width_m)
    if chip_w <= 0.0:
        raise SystemExit("--chip_width_m must be > 0")
    chip_h = chip_w * (dy / dx)  # preserve aspect ratio

    def mapx(x: int) -> float:
        return (x - xmin) / dx * chip_w

    def mapy(y: int) -> float:
        return (y - ymin) / dy * chip_h

    units = sorted(stats.keys())

    # FLP: one rectangle per Bxx, using the bbox of its current source nodes.
    out_flp.parent.mkdir(parents=True, exist_ok=True)
    with open(out_flp, "w") as f:
        f.write("# <unit-name>\t<width>\t<height>\t<left-x>\t<bottom-y>  (meters)\n")
        for b in units:
            bb = stats[b]
            lx = mapx(bb.xmin)
            by = mapy(bb.ymin)
            rx = mapx(bb.xmax)
            ty = mapy(bb.ymax)
            w = max(1e-12, rx - lx)
            h = max(1e-12, ty - by)
            f.write(f"B{b}\t{w:.12g}\t{h:.12g}\t{lx:.12g}\t{by:.12g}\n")

    # PTrace: header (unit names) + one timestep (power in W).
    out_ptrace.parent.mkdir(parents=True, exist_ok=True)
    with open(out_ptrace, "w") as f:
        f.write("\t".join([f"B{b}" for b in units]) + "\t\n")
        powers_w = [args.vdd * stats[b].i_sum_a for b in units]
        f.write("\t".join([f"{p:.12g}" for p in powers_w]) + "\t\n")

    if out_png is not None:
        draw_blocks_png(
            out_png=out_png,
            stats=stats,
            points_by_b=points_by_b,
            mapx=mapx,
            mapy=mapy,
            chip_w=chip_w,
            chip_h=chip_h,
            stride=args.blocks_stride,
        )

    total_i = sum(stats[b].i_sum_a for b in units)
    total_p = args.vdd * total_i
    print(f"Wrote FLP:    {out_flp}")
    print(f"Wrote PTrace: {out_ptrace}")
    if out_png is not None:
        print(f"Wrote blocks PNG: {out_png}")
    print(f"Chip size used in FLP: {chip_w:g} m x {chip_h:g} m")
    print(f"Units: {len(units)}")
    print(f"Total I: {total_i:g} A")
    print(f"Total P: {total_p:g} W (VDD={args.vdd:g} V)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())


