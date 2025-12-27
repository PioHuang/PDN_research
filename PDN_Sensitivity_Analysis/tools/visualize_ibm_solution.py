#!/usr/bin/env python3
"""
Visualize IBM Power Grid Benchmark *.solution (e.g., ibmpg1.solution)

Modes:
  - voltage_heatmap: show absolute voltage for chosen net(s)
  - scatter: scatter plot of absolute voltage
  - combined_ir: binned IR drop % using VDD nets (n1/n3) and GND nets (n0/n2)

Defaults are chosen to mimic VoltSpot's ir_drop visualization:
  cmap=hot, bins=69, vdd=1.8
"""

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np


NODE_RE = re.compile(r"^n(?P<net>\d+)_(?P<x>\d+)_(?P<y>\d+)$")


def parse_solution(path: Path) -> Dict[int, List[Tuple[int, int, float]]]:
    """
    Return net -> list of (x, y, voltage).
    """
    nets: Dict[int, List[Tuple[int, int, float]]] = {}
    with path.open("r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            name = parts[0].strip()
            try:
                v = float(parts[1])
            except ValueError:
                continue
            m = NODE_RE.match(name)
            if not m:
                continue
            net = int(m.group("net"))
            x = int(m.group("x"))
            y = int(m.group("y"))
            nets.setdefault(net, []).append((x, y, v))
    return nets


def _scatter(points: List[Tuple[int, int, float]], title: str, out: Path,
             cmap: str, point_size: float, vmin: float | None, vmax: float | None, no_axes: bool):
    xs = np.array([p[0] for p in points], dtype=float)
    ys = np.array([p[1] for p in points], dtype=float)
    vs = np.array([p[2] for p in points], dtype=float)

    fig, ax = plt.subplots(figsize=(10, 8))
    sc = ax.scatter(xs, ys, c=vs, s=point_size, cmap=cmap, vmin=vmin, vmax=vmax, linewidths=0)
    ax.set_title(title, fontsize=14, fontweight="bold")
    if no_axes:
        ax.set_xticks([]); ax.set_yticks([])
    else:
        ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_aspect("equal", adjustable="box")
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label("Voltage (V)")
    out.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


def _bin_and_heatmap(points: List[Tuple[int, int, float]], title: str, out: Path,
                     bins: int, cmap: str, vmin: float | None, vmax: float | None,
                     no_axes: bool, cbar_label: str):
    xs = np.array([p[0] for p in points], dtype=float)
    ys = np.array([p[1] for p in points], dtype=float)
    vs = np.array([p[2] for p in points], dtype=float)

    x_min, x_max = xs.min(), xs.max()
    y_min, y_max = ys.min(), ys.max()

    sum_grid, _, _ = np.histogram2d(xs, ys, bins=bins, range=[[x_min, x_max], [y_min, y_max]], weights=vs)
    cnt_grid, _, _ = np.histogram2d(xs, ys, bins=bins, range=[[x_min, x_max], [y_min, y_max]])
    with np.errstate(invalid="ignore", divide="ignore"):
        mean_grid = sum_grid / cnt_grid
    mean_grid[cnt_grid == 0] = np.nan

    img = mean_grid.T
    cm = plt.get_cmap(cmap).copy()
    cm.set_bad(color="black")

    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(img, origin="lower", aspect="auto", cmap=cm, vmin=vmin, vmax=vmax)

    finite = img[np.isfinite(img)]
    if finite.size > 0:
        max_val = float(np.max(finite))
        avg_val = float(np.mean(finite))
        ax.set_title(f"{title}\nMax: {max_val:.3f}, Avg: {avg_val:.3f}", fontsize=16, fontweight="bold")
    else:
        ax.set_title(title, fontsize=16, fontweight="bold")

    if no_axes:
        ax.set_xticks([]); ax.set_yticks([])
    else:
        ax.set_xlabel("Column (binned)"); ax.set_ylabel("Row (binned)")

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(cbar_label, fontsize=12)

    out.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


def _combined_ir(nets: Dict[int, List[Tuple[int, int, float]]], out: Path,
                 vdd_ideal: float, bins: int, cmap: str, no_axes: bool):
    # Separate VDD and GND nets
    vdd_pts: List[Tuple[int, int, float]] = []
    gnd_pts: List[Tuple[int, int, float]] = []
    for net_id, pts in nets.items():
        if net_id in (1, 3):
            vdd_pts.extend(pts)
        elif net_id in (0, 2):
            gnd_pts.extend(pts)

    if not vdd_pts or not gnd_pts:
        print("Error: missing VDD or GND points; cannot compute combined IR.")
        return

    v_xs = np.array([p[0] for p in vdd_pts]); v_ys = np.array([p[1] for p in vdd_pts]); v_vs = np.array([p[2] for p in vdd_pts])
    g_xs = np.array([p[0] for p in gnd_pts]); g_ys = np.array([p[1] for p in gnd_pts]); g_vs = np.array([p[2] for p in gnd_pts])

    x_min = min(v_xs.min(), g_xs.min()); x_max = max(v_xs.max(), g_xs.max())
    y_min = min(v_ys.min(), g_ys.min()); y_max = max(v_ys.max(), g_ys.max())

    v_sum, _, _ = np.histogram2d(v_xs, v_ys, bins=bins, range=[[x_min, x_max], [y_min, y_max]], weights=v_vs)
    v_cnt, _, _ = np.histogram2d(v_xs, v_ys, bins=bins, range=[[x_min, x_max], [y_min, y_max]])
    g_sum, _, _ = np.histogram2d(g_xs, g_ys, bins=bins, range=[[x_min, x_max], [y_min, y_max]], weights=g_vs)
    g_cnt, _, _ = np.histogram2d(g_xs, g_ys, bins=bins, range=[[x_min, x_max], [y_min, y_max]])

    def _nan_fill_with_neighbors(arr: np.ndarray, max_iter: int = 100) -> np.ndarray:
        """Fill NaNs by averaging 3x3 neighbors iteratively; fallback to global mean."""
        out_arr = arr.copy()
        for _ in range(max_iter):
            nan_idx = np.argwhere(np.isnan(out_arr))
            if nan_idx.size == 0:
                break
            changed = False
            for i, j in nan_idx:
                r0 = max(0, i - 1); r1 = min(out_arr.shape[0], i + 2)
                c0 = max(0, j - 1); c1 = min(out_arr.shape[1], j + 2)
                window = out_arr[r0:r1, c0:c1]
                finite = window[~np.isnan(window)]
                if finite.size > 0:
                    out_arr[i, j] = float(np.mean(finite))
                    changed = True
            if not changed:
                break
        finite_all = out_arr[~np.isnan(out_arr)]
        if finite_all.size > 0:
            out_arr[np.isnan(out_arr)] = float(np.mean(finite_all))
        return out_arr

    with np.errstate(invalid="ignore", divide="ignore"):
        v_mean = v_sum / v_cnt
        g_mean = g_sum / g_cnt

    v_mean = _nan_fill_with_neighbors(v_mean)
    g_mean = _nan_fill_with_neighbors(g_mean)

    ir = 100.0 * (1.0 - (v_mean - g_mean) / vdd_ideal)

    img = ir.T
    cm = plt.get_cmap(cmap).copy()
    cm.set_bad(color="black")

    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(img, origin="lower", aspect="auto", cmap=cm)

    finite = img[np.isfinite(img)]
    max_val = float(np.max(finite)) if finite.size else 0.0
    avg_val = float(np.mean(finite)) if finite.size else 0.0
    ax.set_title(f"Golden IR Drop (Combined VDD-GND)\nLayer 0\nMax: {max_val:.3f}%, Avg: {avg_val:.3f}%", fontsize=16, fontweight="bold")

    if no_axes:
        ax.set_xticks([]); ax.set_yticks([])
    else:
        ax.set_xlabel("Column (binned)"); ax.set_ylabel("Row (binned)")

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("IR Drop (%)", fontsize=12)

    out.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


def main():
    ap = argparse.ArgumentParser(description="Visualize ibmpg*.solution")
    ap.add_argument("--solution", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--mode", choices=["scatter", "voltage_heatmap", "combined_ir"], default="combined_ir")
    ap.add_argument("--net", type=str, default="all", help="net id or 'all' (for voltage modes)")
    ap.add_argument("--vdd", type=float, default=1.8, help="ideal VDD for IR drop")
    ap.add_argument("--bins", type=int, default=69, help="grid resolution for heatmaps")
    ap.add_argument("--cmap", type=str, default="hot")
    ap.add_argument("--point-size", type=float, default=1.0)
    ap.add_argument("--vmin", type=float, default=None)
    ap.add_argument("--vmax", type=float, default=None)
    ap.add_argument("--no-axes", action="store_true")
    args = ap.parse_args()

    nets = parse_solution(args.solution)
    if not nets:
        print("Error: no data parsed from solution.")
        return

    if args.mode == "combined_ir":
        _combined_ir(nets, args.out, args.vdd, args.bins, args.cmap, args.no_axes)
        return

    # Voltage-based modes
    if args.net != "all":
        try:
            target = int(args.net)
        except ValueError:
            print("Error: --net must be int or 'all'")
            return
        if target not in nets:
            print(f"Error: net {target} not found; available: {sorted(nets.keys())}")
            return
        pts = nets[target]
        label = f"net={target}"
    else:
        pts = [p for pts in nets.values() for p in pts]
        label = "all nets"

    if args.mode == "scatter":
        _scatter(pts, f"{args.solution.name} ({label})", args.out,
                 args.cmap, args.point_size, args.vmin, args.vmax, args.no_axes)
    else:
        _bin_and_heatmap(pts, f"{args.solution.name} ({label})", args.out,
                         args.bins, args.cmap, args.vmin, args.vmax, args.no_axes,
                         "Voltage (V)")


if __name__ == "__main__":
    main()
