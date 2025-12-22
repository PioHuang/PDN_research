#!/usr/bin/env python3
"""
spice_netviz.py

Parse an extracted SPICE PDN netlist (like ibmpg1.spice) and visualize the
per-layer resistor network as an SVG.

Examples:
  python3 PDN_research/PDN_Sensitivity_Analysis/tools/spice_netviz.py --netlist PDN_research/PDN_Sensitivity_Analysis/spice/ibmpg1.spice --list-layers
  python3 PDN_research/PDN_Sensitivity_Analysis/tools/spice_netviz.py --netlist PDN_research/PDN_Sensitivity_Analysis/spice/ibmpg1.spice --layer "M5,VDD" --out m5_vdd.svg
  python3 PDN_research/PDN_Sensitivity_Analysis/tools/spice_netviz.py --netlist PDN_research/PDN_Sensitivity_Analysis/spice/ibmpg1.spice --all --outdir out_svg
  python3 PDN_research/PDN_Sensitivity_Analysis/tools/spice_netviz.py --netlist PDN_research/PDN_Sensitivity_Analysis/spice/ibmpg1.spice --layer "M6,GND" --bbox 0 10000 0 10000 --out crop.svg
  python3 PDN_research/PDN_Sensitivity_Analysis/tools/spice_netviz.py --netlist PDN_research/PDN_Sensitivity_Analysis/spice/ibmpg1.spice --out-3d-png pdn_3d.png
  python3 PDN_research/PDN_Sensitivity_Analysis/tools/spice_netviz.py --netlist PDN_research/PDN_Sensitivity_Analysis/spice/ibmpg1.spice --show-3d --3d-stride 5
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Optional


LAYER_LINE_RE = re.compile(r"^\*\s*layer:\s*(?P<label>.+?)\s*$", re.IGNORECASE)
NET_ID_RE = re.compile(r"\bnet:\s*(?P<net>\d+)\b", re.IGNORECASE)
METAL_RE = re.compile(r"\bM(?P<num>\d+)\b", re.IGNORECASE)

# Matches "n1_333_383" and "_X_n1_333_383" forms.
NODE_WITH_XY_RE = re.compile(r"^(?:_X_)?(?P<prefix>n\d+)_(?P<x>\d+)_(?P<y>\d+)$")


@dataclass(frozen=True)
class Edge:
    n1: str
    n2: str
    x1: int
    y1: int
    x2: int
    y2: int
    r_ohm: float


@dataclass(frozen=True)
class LayerSummary:
    label: str
    start_line: int
    resistor_edges: int


@dataclass(frozen=True)
class Segment3D:
    x1: float
    y1: float
    z1: float
    x2: float
    y2: float
    z2: float
    group: str


@dataclass(frozen=True)
class Branch:
    layer: int
    row: int
    col: int
    direction: str  # "x" or "y"
    resistance: float
    node1: int
    node2: int


def _parse_xy(node: str) -> Optional[tuple[int, int]]:
    match = NODE_WITH_XY_RE.match(node)
    if not match:
        return None
    return int(match.group("x")), int(match.group("y"))


def _is_resistor_line(tokens: list[str]) -> bool:
    return bool(tokens) and tokens[0] and tokens[0][0] in {"R", "r"} and len(tokens) >= 4


def _safe_float(text: str) -> Optional[float]:
    try:
        return float(text)
    except ValueError:
        return None


def iter_netlist_lines(path: Path) -> Iterator[tuple[int, str]]:
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for lineno, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            yield lineno, line


def scan_layers(path: Path) -> list[LayerSummary]:
    current_label: Optional[str] = None
    layer_start: int = 0
    layer_counts: dict[str, int] = {}
    layer_start_line: dict[str, int] = {}

    for lineno, line in iter_netlist_lines(path):
        match = LAYER_LINE_RE.match(line)
        if match:
            current_label = match.group("label")
            layer_start = lineno
            layer_start_line.setdefault(current_label, layer_start)
            layer_counts.setdefault(current_label, 0)
            continue

        if current_label is None:
            continue

        tokens = line.split()
        if _is_resistor_line(tokens):
            n1, n2 = tokens[1], tokens[2]
            if n1.startswith("n") and n2.startswith("n"):
                layer_counts[current_label] += 1

    summaries: list[LayerSummary] = []
    for label, count in layer_counts.items():
        summaries.append(
            LayerSummary(label=label, start_line=layer_start_line.get(label, 0), resistor_edges=count)
        )

    summaries.sort(key=lambda s: s.start_line)
    return summaries


def _matches_layer(label: str, wanted: str) -> bool:
    if label == wanted:
        return True
    return wanted.lower() in label.lower()


def _extract_net_id(label: str) -> Optional[int]:
    match = NET_ID_RE.search(label)
    if not match:
        return None
    return int(match.group("net"))


def _extract_metal_num(label: str) -> Optional[int]:
    match = METAL_RE.search(label)
    if not match:
        return None
    return int(match.group("num"))


def _net_id_from_node(node: str) -> Optional[int]:
    if not node.startswith("n"):
        return None
    underscore = node.find("_")
    if underscore <= 1:
        return None
    try:
        return int(node[1:underscore])
    except ValueError:
        return None


def build_net_to_metal_map(path: Path) -> dict[int, int]:
    net_to_metal: dict[int, int] = {}
    for _lineno, line in iter_netlist_lines(path):
        match = LAYER_LINE_RE.match(line)
        if not match:
            continue
        label = match.group("label")
        net_id = _extract_net_id(label)
        metal = _extract_metal_num(label)
        if net_id is None or metal is None:
            continue
        net_to_metal[net_id] = metal
    return net_to_metal


def extract_layer_edges(
    path: Path,
    *,
    layer_wanted: str,
    bbox: Optional[tuple[int, int, int, int]],
    min_r: Optional[float],
    max_r: Optional[float],
    max_edges: Optional[int],
) -> tuple[str, list[Edge]]:
    current_label: Optional[str] = None
    selected_label: Optional[str] = None
    edges: list[Edge] = []

    for _lineno, line in iter_netlist_lines(path):
        match = LAYER_LINE_RE.match(line)
        if match:
            current_label = match.group("label")
            continue

        if current_label is None:
            continue

        if not _matches_layer(current_label, layer_wanted):
            continue

        selected_label = current_label

        tokens = line.split()
        if not _is_resistor_line(tokens):
            continue

        n1, n2 = tokens[1], tokens[2]
        if not (n1.startswith("n") and n2.startswith("n")):
            continue

        r_ohm = _safe_float(tokens[3])
        if r_ohm is None:
            continue

        if min_r is not None and r_ohm < min_r:
            continue
        if max_r is not None and r_ohm > max_r:
            continue

        xy1 = _parse_xy(n1)
        xy2 = _parse_xy(n2)
        if xy1 is None or xy2 is None:
            continue
        x1, y1 = xy1
        x2, y2 = xy2

        if bbox is not None:
            xmin, xmax, ymin, ymax = bbox
            if not (xmin <= x1 <= xmax and xmin <= x2 <= xmax and ymin <= y1 <= ymax and ymin <= y2 <= ymax):
                continue

        edges.append(Edge(n1=n1, n2=n2, x1=x1, y1=y1, x2=x2, y2=y2, r_ohm=r_ohm))
        if max_edges is not None and len(edges) >= max_edges:
            break

    if selected_label is None:
        available = ", ".join(s.label for s in scan_layers(path))
        raise SystemExit(
            f"No matching layer for {layer_wanted!r}. Available layers: {available or '(none found)'}"
        )

    return selected_label, edges


def extract_3d_segments(
    path: Path,
    *,
    layer_selector: Optional[str],
    bbox: Optional[tuple[int, int, int, int]],
    min_r: Optional[float],
    max_r: Optional[float],
    max_edges: Optional[int],
    include_vias: bool,
    z_step: float,
) -> list[Segment3D]:
    net_to_metal = build_net_to_metal_map(path)
    if not net_to_metal:
        raise SystemExit("No parseable '* layer:' headers (need 'M#' and 'net: <id>') for 3D export.")

    metal_min = min(net_to_metal.values())

    def z_for_net(net_id: int) -> float:
        metal = net_to_metal.get(net_id)
        if metal is None:
            return 0.0
        return float(metal - metal_min) * z_step

    selected_nets: Optional[set[int]] = None
    if layer_selector is not None:
        selected_nets = set()
        for _lineno, line in iter_netlist_lines(path):
            match = LAYER_LINE_RE.match(line)
            if not match:
                continue
            label = match.group("label")
            if not _matches_layer(label, layer_selector):
                continue
            net_id = _extract_net_id(label)
            if net_id is not None:
                selected_nets.add(net_id)

    segments: list[Segment3D] = []
    current_label: Optional[str] = None
    current_net: Optional[int] = None

    for _lineno, line in iter_netlist_lines(path):
        match = LAYER_LINE_RE.match(line)
        if match:
            current_label = match.group("label")
            current_net = _extract_net_id(current_label)
            continue

        if selected_nets is not None:
            if current_net is None or current_net not in selected_nets:
                continue

        tokens = line.split()
        if not _is_resistor_line(tokens):
            continue

        n1, n2 = tokens[1], tokens[2]
        r_ohm = _safe_float(tokens[3])
        if r_ohm is None:
            continue
        if min_r is not None and r_ohm < min_r:
            continue
        if max_r is not None and r_ohm > max_r:
            continue

        xy1 = _parse_xy(n1)
        xy2 = _parse_xy(n2)
        if xy1 is None or xy2 is None:
            continue
        x1, y1 = xy1
        x2, y2 = xy2
        if bbox is not None:
            xmin, xmax, ymin, ymax = bbox
            if not (xmin <= x1 <= xmax and xmin <= x2 <= xmax and ymin <= y1 <= ymax and ymin <= y2 <= ymax):
                continue

        net1 = _net_id_from_node(n1)
        net2 = _net_id_from_node(n2)
        if net1 is None or net2 is None or net1 != net2:
            continue
        if net1 not in net_to_metal:
            continue

        group = current_label or f"net:{net1}"
        z = z_for_net(net1)
        segments.append(Segment3D(x1=float(x1), y1=float(y1), z1=z, x2=float(x2), y2=float(y2), z2=z, group=group))
        if max_edges is not None and len(segments) >= max_edges:
            break

    if include_vias:
        for _lineno, line in iter_netlist_lines(path):
            tokens = line.split()
            if not tokens or tokens[0][0] not in {"V", "v"}:
                continue
            if len(tokens) < 4:
                continue
            a, b = tokens[1], tokens[2]
            v = _safe_float(tokens[3])
            if v is None or abs(v) > 1e-15:
                continue

            # Via-like sources: connect two grid nodes at same (x,y).
            if not (a.startswith("n") and b.startswith("n")):
                continue
            xy_a = _parse_xy(a)
            xy_b = _parse_xy(b)
            if xy_a is None or xy_b is None:
                continue
            xa, ya = xy_a
            xb, yb = xy_b
            if xa != xb or ya != yb:
                continue

            if bbox is not None:
                xmin, xmax, ymin, ymax = bbox
                if not (xmin <= xa <= xmax and ymin <= ya <= ymax):
                    continue

            net_a = _net_id_from_node(a)
            net_b = _net_id_from_node(b)
            if net_a is None or net_b is None:
                continue
            if net_a not in net_to_metal or net_b not in net_to_metal:
                continue
            if selected_nets is not None and (net_a not in selected_nets and net_b not in selected_nets):
                continue

            z1 = z_for_net(net_a)
            z2 = z_for_net(net_b)
            if z1 == z2:
                continue
            segments.append(Segment3D(x1=float(xa), y1=float(ya), z1=z1, x2=float(xa), y2=float(ya), z2=z2, group="vias"))

    return segments


def _ensure_mpl_cache_dir() -> None:
    if os.environ.get("MPLCONFIGDIR"):
        return
    base = Path(os.environ.get("TMPDIR", "/tmp"))
    cache_dir = base / "mplconfig"
    try:
        cache_dir.mkdir(parents=True, exist_ok=True)
    except OSError:
        return
    os.environ["MPLCONFIGDIR"] = str(cache_dir)


def read_branches_csv(path: Path) -> list[Branch]:
    if not path.exists():
        raise SystemExit(f"branches.csv not found: {path}")

    branches: list[Branch] = []
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"layer", "row", "col", "direction", "resistance", "node1", "node2"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise SystemExit(f"Invalid branches.csv header; expected at least {sorted(required)}")

        for rec in reader:
            try:
                direction = str(rec["direction"]).strip().lower()
                if direction not in {"x", "y"}:
                    continue
                branches.append(
                    Branch(
                        layer=int(rec["layer"]),
                        row=int(rec["row"]),
                        col=int(rec["col"]),
                        direction=direction,
                        resistance=float(rec["resistance"]),
                        node1=int(rec["node1"]),
                        node2=int(rec["node2"]),
                    )
                )
            except (KeyError, ValueError, TypeError):
                continue

    if not branches:
        raise SystemExit(f"No branches parsed from {path}")
    return branches


def plot_branches_virtual_grid(
    *,
    branches_csv: Path,
    out_png: Optional[Path],
    show: bool,
    stride: int,
    lw: float,
    alpha: float,
    dpi: int,
) -> None:
    if stride < 1:
        raise SystemExit("--branches-stride must be >= 1")

    branches = read_branches_csv(branches_csv)

    max_row = 0
    max_col = 0
    x_segments: list[tuple[tuple[float, float], tuple[float, float]]] = []
    y_segments: list[tuple[tuple[float, float], tuple[float, float]]] = []

    for b in branches:
        if b.direction == "x":
            x_segments.append(((float(b.col), float(b.row)), (float(b.col + 1), float(b.row))))
            max_col = max(max_col, b.col + 1)
            max_row = max(max_row, b.row)
        else:
            y_segments.append(((float(b.col), float(b.row)), (float(b.col), float(b.row + 1))))
            max_col = max(max_col, b.col)
            max_row = max(max_row, b.row + 1)

    _ensure_mpl_cache_dir()

    import matplotlib

    if not show:
        matplotlib.use("Agg")

    from matplotlib import pyplot as plt  # noqa: E402
    from matplotlib.collections import LineCollection  # noqa: E402

    fig, ax = plt.subplots(figsize=(9.0, 9.0), dpi=dpi)
    ax.set_aspect("equal", adjustable="box")

    x_plot = x_segments[::stride]
    y_plot = y_segments[::stride]

    if x_plot:
        ax.add_collection(LineCollection(x_plot, colors="tab:red", linewidths=lw, alpha=alpha))
    if y_plot:
        ax.add_collection(LineCollection(y_plot, colors="tab:blue", linewidths=lw, alpha=alpha))

    ax.set_xlim(-0.5, float(max_col) + 0.5)
    ax.set_ylim(-0.5, float(max_row) + 0.5)
    ax.set_xlabel("col")
    ax.set_ylabel("row")
    ax.set_title(
        f"Virtual grid branches ({branches_csv.name})  "
        f"edges: x={len(x_segments)} y={len(y_segments)}  "
        f"plotted: x={len(x_plot)} y={len(y_plot)}"
    )
    ax.grid(False)

    fig.tight_layout()
    if out_png is not None:
        out_png.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_png)

    if show:
        plt.show()

    plt.close(fig)


def render_3d_matplotlib(
    *,
    segments: list[Segment3D],
    elev: float,
    azim: float,
    stride: int,
    linewidth: float,
    alpha: float,
    vias_stride: int,
    vias_linewidth: float,
    vias_alpha: float,
    figsize: tuple[float, float],
    dpi: int,
    interactive: bool,
):
    if not segments:
        raise SystemExit("No 3D segments to plot (check filters).")
    if stride < 1:
        raise SystemExit("--3d-stride must be >= 1")
    if vias_stride < 1:
        raise SystemExit("--vias-stride must be >= 1")

    _ensure_mpl_cache_dir()

    import matplotlib

    if not interactive:
        matplotlib.use("Agg")

    from matplotlib import pyplot as plt  # noqa: E402
    from mpl_toolkits.mplot3d.art3d import Line3DCollection  # noqa: E402

    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111, projection="3d")
    ax.set_axis_off()
    ax.view_init(elev=elev, azim=azim)

    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []
    for s in segments:
        xs.extend([s.x1, s.x2])
        ys.extend([s.y1, s.y2])
        zs.extend([s.z1, s.z2])

    ax.set_xlim(min(xs), max(xs))
    ax.set_ylim(min(ys), max(ys))
    ax.set_zlim(min(zs), max(zs))

    # Vias (vertical lines).
    via_lines = [((s.x1, s.y1, s.z1), (s.x2, s.y2, s.z2)) for s in segments if s.group == "vias"]
    if via_lines:
        via_lines = via_lines[::vias_stride]
        lc_vias = Line3DCollection(via_lines, colors="#111111", linewidths=vias_linewidth, alpha=vias_alpha)
        ax.add_collection3d(lc_vias)

    # Metal lines colored by z-plane.
    metal_segments = [s for s in segments if s.group != "vias"]
    metal_segments = metal_segments[::stride]
    if metal_segments:
        z_values = sorted({s.z1 for s in metal_segments})
        cmap = plt.get_cmap("tab10")
        z_to_color = {z: cmap(i % cmap.N) for i, z in enumerate(z_values)}
        metal_lines = [((s.x1, s.y1, s.z1), (s.x2, s.y2, s.z2)) for s in metal_segments]
        metal_colors = [z_to_color.get(s.z1, "#1f77b4") for s in metal_segments]
        lc_metal = Line3DCollection(metal_lines, colors=metal_colors, linewidths=linewidth, alpha=alpha)
        ax.add_collection3d(lc_metal)

    return fig, ax


def write_3d_png(
    out_path: Path,
    *,
    segments: list[Segment3D],
    elev: float,
    azim: float,
    stride: int,
    linewidth: float,
    alpha: float,
    vias_stride: int,
    vias_linewidth: float,
    vias_alpha: float,
    dpi: int,
    figsize: tuple[float, float],
) -> None:
    fig, _ax = render_3d_matplotlib(
        segments=segments,
        elev=elev,
        azim=azim,
        stride=stride,
        linewidth=linewidth,
        alpha=alpha,
        vias_stride=vias_stride,
        vias_linewidth=vias_linewidth,
        vias_alpha=vias_alpha,
        figsize=figsize,
        dpi=dpi,
        interactive=False,
    )
    from matplotlib import pyplot as plt  # noqa: E402

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(pad=0)
    fig.savefig(out_path, bbox_inches="tight", pad_inches=0)
    plt.close(fig)


def _rgb_to_hex(rgb: tuple[int, int, int]) -> str:
    r, g, b = rgb
    return f"#{r:02x}{g:02x}{b:02x}"


def _lerp(a: float, b: float, t: float) -> float:
    return a + (b - a) * t


def _lerp_rgb(c1: tuple[int, int, int], c2: tuple[int, int, int], t: float) -> tuple[int, int, int]:
    t = max(0.0, min(1.0, t))
    return (
        int(round(_lerp(c1[0], c2[0], t))),
        int(round(_lerp(c1[1], c2[1], t))),
        int(round(_lerp(c1[2], c2[2], t))),
    )


def _color_for_r(r_ohm: float, r_min: float, r_max: float) -> str:
    if r_ohm <= 0:
        return "#000000"
    a = math.log10(r_min) if r_min > 0 else math.log10(min(v for v in [r_max, 1e-12] if v > 0))
    b = math.log10(r_max) if r_max > 0 else a + 1.0
    x = math.log10(r_ohm)
    t = 0.0 if b == a else (x - a) / (b - a)
    t = max(0.0, min(1.0, t))

    # blue -> yellow -> red
    blue = (44, 123, 182)
    yellow = (255, 255, 191)
    red = (215, 25, 28)
    if t < 0.5:
        return _rgb_to_hex(_lerp_rgb(blue, yellow, t / 0.5))
    return _rgb_to_hex(_lerp_rgb(yellow, red, (t - 0.5) / 0.5))


def write_svg(
    out_path: Path,
    *,
    title: str,
    edges: list[Edge],
    width_px: int,
    height_px: int,
    margin_px: int,
) -> None:
    if not edges:
        raise SystemExit("No edges to plot after filtering.")

    xs = [e.x1 for e in edges] + [e.x2 for e in edges]
    ys = [e.y1 for e in edges] + [e.y2 for e in edges]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)

    dx = max(1, xmax - xmin)
    dy = max(1, ymax - ymin)

    sx = (width_px - 2 * margin_px) / dx
    sy = (height_px - 2 * margin_px) / dy
    scale = min(sx, sy)

    def to_px(x: int, y: int) -> tuple[float, float]:
        px = margin_px + (x - xmin) * scale
        py = height_px - (margin_px + (y - ymin) * scale)  # flip y for usual screen coords
        return px, py

    r_min = min(e.r_ohm for e in edges if e.r_ohm > 0)
    r_max = max(e.r_ohm for e in edges if e.r_ohm > 0)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as handle:
        handle.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        handle.write(
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width_px}" height="{height_px}" '
            f'viewBox="0 0 {width_px} {height_px}">\n'
        )
        handle.write('<rect x="0" y="0" width="100%" height="100%" fill="white"/>\n')
        handle.write(
            f'<text x="{margin_px}" y="{margin_px}" font-size="14" '
            f'font-family="monospace" fill="#111">{_escape_xml(title)}</text>\n'
        )
        handle.write('<g stroke-linecap="round" stroke-opacity="0.85">\n')

        # Constant stroke width keeps files smaller and renders faster.
        stroke_width = 1.0
        for edge in edges:
            x1, y1 = to_px(edge.x1, edge.y1)
            x2, y2 = to_px(edge.x2, edge.y2)
            color = _color_for_r(edge.r_ohm, r_min, r_max)
            handle.write(
                f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" '
                f'stroke="{color}" stroke-width="{stroke_width}"/>\n'
            )
        handle.write("</g>\n</svg>\n")


def _escape_xml(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&apos;")
    )


def export_csv(out_dir: Path, *, label: str, edges: list[Edge]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    edges_path = out_dir / f"{_slug(label)}.edges.csv"
    with edges_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["n1", "n2", "x1", "y1", "x2", "y2", "r_ohm"])
        for e in edges:
            writer.writerow([e.n1, e.n2, e.x1, e.y1, e.x2, e.y2, f"{e.r_ohm:.12g}"])


def _slug(text: str) -> str:
    safe = re.sub(r"[^a-zA-Z0-9]+", "_", text).strip("_")
    return safe or "layer"


def _parse_bbox(args: list[str]) -> tuple[int, int, int, int]:
    if len(args) != 4:
        raise SystemExit("--bbox requires 4 integers: xmin xmax ymin ymax")
    xmin, xmax, ymin, ymax = (int(v) for v in args)
    if xmin > xmax or ymin > ymax:
        raise SystemExit("--bbox must satisfy xmin<=xmax and ymin<=ymax")
    return xmin, xmax, ymin, ymax


def segment_stats(segments: list[Segment3D]) -> dict[str, int]:
    stats: dict[str, int] = {"total": len(segments), "vias": 0, "metal": 0}
    for s in segments:
        if s.group == "vias":
            stats["vias"] += 1
        else:
            stats["metal"] += 1
    return stats


def choose_pitch_from_bbox(xmin: float, xmax: float, ymin: float, ymax: float, *, target_pixels: int) -> float:
    span = max(xmax - xmin, ymax - ymin)
    if span <= 0:
        return 1.0
    raw = span / float(max(1, target_pixels))
    # Round to a stable "nice" value without overshooting too much.
    # With typical extracted PDNs, overshooting makes many edges span multiple pixels and hurts accuracy.
    step = 25.0 if raw >= 100 else 10.0
    return max(step, round(raw / step) * step)


def _iter_nodes_from_edges(edges: list[Edge]) -> Iterator[tuple[int, int]]:
    for e in edges:
        yield e.x1, e.y1
        yield e.x2, e.y2


def _bbox_from_edges(edges: list[Edge]) -> tuple[int, int, int, int]:
    xs: list[int] = []
    ys: list[int] = []
    for x, y in _iter_nodes_from_edges(edges):
        xs.append(x)
        ys.append(y)
    if not xs:
        return 0, 0, 0, 0
    return min(xs), max(xs), min(ys), max(ys)


def _pixel_index(x: int, y: int, *, xmin: int, ymin: int, pitch: float) -> tuple[int, int]:
    px = int(math.floor((x - xmin) / pitch))
    py = int(math.floor((y - ymin) / pitch))
    return px, py


def _node_id(px: int, py: int, nx: int) -> int:
    return py * nx + px


def _cg_solve(
    matvec,
    b: list[float],
    *,
    x0: Optional[list[float]] = None,
    max_iter: int = 2000,
    tol: float = 1e-8,
) -> list[float]:
    n = len(b)
    x = [0.0] * n if x0 is None else list(x0)
    r = [bi - ai for bi, ai in zip(b, matvec(x))]
    p = list(r)

    def dot(u: list[float], v: list[float]) -> float:
        return sum(ui * vi for ui, vi in zip(u, v))

    rsold = dot(r, r)
    if rsold == 0.0:
        return x

    bnorm = math.sqrt(dot(b, b))
    if bnorm == 0.0:
        return x

    for _ in range(max_iter):
        Ap = matvec(p)
        denom = dot(p, Ap)
        if denom == 0.0:
            break
        alpha = rsold / denom
        x = [xi + alpha * pi for xi, pi in zip(x, p)]
        r = [ri - alpha * api for ri, api in zip(r, Ap)]
        rsnew = dot(r, r)
        if math.sqrt(rsnew) / bnorm < tol:
            break
        beta = rsnew / rsold
        p = [ri + beta * pi for ri, pi in zip(r, p)]
        rsold = rsnew

    return x


def _estimate_pixel_axis_r(
    *,
    nx: int,
    ny: int,
    links: dict[tuple[int, int], float],
    axis: str,
    ridge: float = 1e-12,
) -> list[float]:
    """
    Solve Rx or Ry from neighbor link resistances using weighted least squares:
      0.5*R[p] + 0.5*R[q] ~= Rlink(p,q)
    where weights are link conductances.
    """
    n = nx * ny
    edges: list[tuple[int, int, float, float]] = []
    for (a, b), g in links.items():
        if g <= 0:
            continue
        rlink = 1.0 / g
        edges.append((a, b, g, rlink))

    if not edges:
        return [float("nan")] * n

    acoef = 0.5

    def matvec(x: list[float]) -> list[float]:
        y = [ridge * xi for xi in x]
        for i, j, w, _rlink in edges:
            t = w * (acoef * x[i] + acoef * x[j])
            y[i] += acoef * t
            y[j] += acoef * t
        return y

    rhs = [0.0] * n
    for i, j, w, rlink in edges:
        t = w * rlink
        rhs[i] += acoef * t
        rhs[j] += acoef * t

    x = _cg_solve(matvec, rhs, max_iter=2000, tol=1e-9)
    x = [xi if xi > 0 else 1e-18 for xi in x]
    return x


def export_pixel_model_vdd(
    netlist: Path,
    out_dir: Path,
    *,
    pitch: Optional[float],
    target_pixels: int,
    rz_per_via: float,
    bbox: Optional[tuple[int, int, int, int]],
    max_edges: Optional[int],
) -> None:
    """
    Export a simple 2-layer (M5/M6) VDD pixel model:
      - per-layer Rx,Ry (per pixel)
      - per-pixel via count and Rz (between M5 and M6 at that pixel)
    """
    label_m5, edges_m5 = extract_layer_edges(
        netlist, layer_wanted="M5,VDD", bbox=bbox, min_r=None, max_r=None, max_edges=max_edges
    )
    label_m6, edges_m6 = extract_layer_edges(
        netlist, layer_wanted="M6,VDD", bbox=bbox, min_r=None, max_r=None, max_edges=max_edges
    )

    xmin1, xmax1, ymin1, ymax1 = _bbox_from_edges(edges_m5)
    xmin2, xmax2, ymin2, ymax2 = _bbox_from_edges(edges_m6)
    xmin = min(xmin1, xmin2)
    ymin = min(ymin1, ymin2)
    xmax = max(xmax1, xmax2)
    ymax = max(ymax1, ymax2)

    if pitch is None:
        pitch = choose_pitch_from_bbox(xmin, xmax, ymin, ymax, target_pixels=target_pixels)

    if pitch <= 0:
        raise SystemExit("--pitch must be > 0")

    nx = int(math.floor((xmax - xmin) / pitch)) + 1
    ny = int(math.floor((ymax - ymin) / pitch)) + 1
    if nx <= 0 or ny <= 0:
        raise SystemExit("Invalid grid size from bbox/pitch.")

    def accumulate_links(edges: list[Edge]) -> tuple[dict[tuple[int, int], float], dict[tuple[int, int], float], dict[str, int]]:
        links_x: dict[tuple[int, int], float] = {}
        links_y: dict[tuple[int, int], float] = {}
        stats = {"kept_neighbor": 0, "skipped_same_pixel": 0, "split_non_neighbor": 0, "skipped_oob": 0}

        for e in edges:
            p1 = _pixel_index(e.x1, e.y1, xmin=xmin, ymin=ymin, pitch=pitch)
            p2 = _pixel_index(e.x2, e.y2, xmin=xmin, ymin=ymin, pitch=pitch)
            if p1 == p2:
                stats["skipped_same_pixel"] += 1
                continue
            x1, y1 = p1
            x2, y2 = p2
            if not (0 <= x1 < nx and 0 <= y1 < ny and 0 <= x2 < nx and 0 <= y2 < ny):
                stats["skipped_oob"] += 1
                continue
            dx = x2 - x1
            dy = y2 - y1
            if e.r_ohm <= 0:
                continue
            if abs(dx) + abs(dy) == 1:
                g = 1.0 / e.r_ohm
                a = _node_id(x1, y1, nx)
                b = _node_id(x2, y2, nx)
                if a > b:
                    a, b = b, a
                    dx, dy = -dx, -dy
                key = (a, b)
                if dy == 0:
                    links_x[key] = links_x.get(key, 0.0) + g
                else:
                    links_y[key] = links_y.get(key, 0.0) + g
                stats["kept_neighbor"] += 1
                continue

            # Split long edges across the pixel grid so the coarse model remains a neighbor mesh.
            steps = abs(dx) + abs(dy)
            if steps <= 0:
                continue
            r_step = e.r_ohm / float(steps)
            g_step = 1.0 / r_step
            cx, cy = x1, y1
            sx = 1 if dx > 0 else (-1 if dx < 0 else 0)
            sy = 1 if dy > 0 else (-1 if dy < 0 else 0)
            for _ in range(abs(dx)):
                nxp, nyp = cx + sx, cy
                a = _node_id(cx, cy, nx)
                b = _node_id(nxp, nyp, nx)
                if a > b:
                    a, b = b, a
                links_x[(a, b)] = links_x.get((a, b), 0.0) + g_step
                cx = nxp
            for _ in range(abs(dy)):
                nxp, nyp = cx, cy + sy
                a = _node_id(cx, cy, nx)
                b = _node_id(nxp, nyp, nx)
                if a > b:
                    a, b = b, a
                links_y[(a, b)] = links_y.get((a, b), 0.0) + g_step
                cy = nyp
            stats["split_non_neighbor"] += 1
        return links_x, links_y, stats

    links_m5_x, links_m5_y, stats_m5 = accumulate_links(edges_m5)
    links_m6_x, links_m6_y, stats_m6 = accumulate_links(edges_m6)

    rx_m5 = _estimate_pixel_axis_r(nx=nx, ny=ny, links=links_m5_x, axis="x")
    ry_m5 = _estimate_pixel_axis_r(nx=nx, ny=ny, links=links_m5_y, axis="y")
    rx_m6 = _estimate_pixel_axis_r(nx=nx, ny=ny, links=links_m6_x, axis="x")
    ry_m6 = _estimate_pixel_axis_r(nx=nx, ny=ny, links=links_m6_y, axis="y")

    # Via counts for VDD: between net 1 and net 3 at same (x,y) using 0V sources.
    via_count = [0] * (nx * ny)
    for _lineno, line in iter_netlist_lines(netlist):
        tokens = line.split()
        if not tokens or tokens[0][0] not in {"V", "v"}:
            continue
        if len(tokens) < 4:
            continue
        a, b = tokens[1], tokens[2]
        v = _safe_float(tokens[3])
        if v is None or abs(v) > 1e-15:
            continue
        if not (a.startswith("n") and b.startswith("n")):
            continue
        net_a = _net_id_from_node(a)
        net_b = _net_id_from_node(b)
        if net_a is None or net_b is None:
            continue
        if {net_a, net_b} != {1, 3}:
            continue
        xy_a = _parse_xy(a)
        xy_b = _parse_xy(b)
        if xy_a is None or xy_b is None:
            continue
        xa, ya = xy_a
        xb, yb = xy_b
        if xa != xb or ya != yb:
            continue
        if bbox is not None:
            bxmin, bxmax, bymin, bymax = bbox
            if not (bxmin <= xa <= bxmax and bymin <= ya <= bymax):
                continue
        px, py = _pixel_index(xa, ya, xmin=xmin, ymin=ymin, pitch=pitch)
        if not (0 <= px < nx and 0 <= py < ny):
            continue
        via_count[_node_id(px, py, nx)] += 1

    rz = [float("inf")] * (nx * ny)
    if rz_per_via <= 0:
        raise SystemExit("--rz-per-via must be > 0")
    for idx, c in enumerate(via_count):
        if c > 0:
            rz[idx] = rz_per_via / float(c)

    out_dir.mkdir(parents=True, exist_ok=True)
    meta = {
        "net": "VDD",
        "layers": {"m5": label_m5, "m6": label_m6},
        "pitch": pitch,
        "origin": {"xmin": xmin, "ymin": ymin},
        "bbox": {"xmin": xmin, "xmax": xmax, "ymin": ymin, "ymax": ymax},
        "grid": {"nx": nx, "ny": ny},
        "rz_per_via": rz_per_via,
        "edge_stats": {"m5": stats_m5, "m6": stats_m6},
    }
    (out_dir / "meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    def write_matrix(path: Path, data: list[float]) -> None:
        with path.open("w", encoding="utf-8", newline="") as handle:
            w = csv.writer(handle)
            for py in range(ny):
                row = data[py * nx : (py + 1) * nx]
                w.writerow(row)

    def write_int_matrix(path: Path, data: list[int]) -> None:
        with path.open("w", encoding="utf-8", newline="") as handle:
            w = csv.writer(handle)
            for py in range(ny):
                row = data[py * nx : (py + 1) * nx]
                w.writerow(row)

    write_matrix(out_dir / "rx_m5.csv", rx_m5)
    write_matrix(out_dir / "ry_m5.csv", ry_m5)
    write_matrix(out_dir / "rx_m6.csv", rx_m6)
    write_matrix(out_dir / "ry_m6.csv", ry_m6)
    write_int_matrix(out_dir / "via_count.csv", via_count)
    write_matrix(out_dir / "rz_via.csv", rz)


def export_pixel_model_vdd_merged(
    netlist: Path,
    out_dir: Path,
    *,
    pitch: Optional[float],
    target_pixels: int,
    rz_per_via: float,
    bbox: Optional[tuple[int, int, int, int]],
    max_edges: Optional[int],
) -> None:
    """
    Export a merged single-layer VDD pixel model.

    Mapping:
      - Build neighbor link conductances on M5 and M6 (splitting long edges across pixels).
      - Merge in-plane conductances by parallel sum: G_total = G_m5 + G_m6.
      - Solve per-pixel Rx/Ry from neighbor links: 0.5*R[p] + 0.5*R[q] ~= Rlink(p,q).
      - Use via_count between net 1 and net 3 to build Rz = rz_per_via / via_count (optional vertical coupling).

    This is a good approximation when M5/M6 are strongly tied by many vias.
    """
    label_m5, edges_m5 = extract_layer_edges(
        netlist, layer_wanted="M5,VDD", bbox=bbox, min_r=None, max_r=None, max_edges=max_edges
    )
    label_m6, edges_m6 = extract_layer_edges(
        netlist, layer_wanted="M6,VDD", bbox=bbox, min_r=None, max_r=None, max_edges=max_edges
    )

    xmin1, xmax1, ymin1, ymax1 = _bbox_from_edges(edges_m5)
    xmin2, xmax2, ymin2, ymax2 = _bbox_from_edges(edges_m6)
    xmin = min(xmin1, xmin2)
    ymin = min(ymin1, ymin2)
    xmax = max(xmax1, xmax2)
    ymax = max(ymax1, ymax2)

    if pitch is None:
        pitch = choose_pitch_from_bbox(xmin, xmax, ymin, ymax, target_pixels=target_pixels)
    if pitch <= 0:
        raise SystemExit("--pitch must be > 0")

    nx = int(math.floor((xmax - xmin) / pitch)) + 1
    ny = int(math.floor((ymax - ymin) / pitch)) + 1
    if nx <= 0 or ny <= 0:
        raise SystemExit("Invalid grid size from bbox/pitch.")

    def accumulate_links(edges: list[Edge]) -> tuple[dict[tuple[int, int], float], dict[tuple[int, int], float], dict[str, int]]:
        links_x: dict[tuple[int, int], float] = {}
        links_y: dict[tuple[int, int], float] = {}
        stats = {"kept_neighbor": 0, "skipped_same_pixel": 0, "split_non_neighbor": 0, "skipped_oob": 0}

        for e in edges:
            p1 = _pixel_index(e.x1, e.y1, xmin=xmin, ymin=ymin, pitch=pitch)
            p2 = _pixel_index(e.x2, e.y2, xmin=xmin, ymin=ymin, pitch=pitch)
            if p1 == p2:
                stats["skipped_same_pixel"] += 1
                continue
            x1, y1 = p1
            x2, y2 = p2
            if not (0 <= x1 < nx and 0 <= y1 < ny and 0 <= x2 < nx and 0 <= y2 < ny):
                stats["skipped_oob"] += 1
                continue
            dx = x2 - x1
            dy = y2 - y1
            if e.r_ohm <= 0:
                continue

            if abs(dx) + abs(dy) == 1:
                g = 1.0 / e.r_ohm
                a = _node_id(x1, y1, nx)
                b = _node_id(x2, y2, nx)
                if a > b:
                    a, b = b, a
                    dx, dy = -dx, -dy
                key = (a, b)
                if dy == 0:
                    links_x[key] = links_x.get(key, 0.0) + g
                else:
                    links_y[key] = links_y.get(key, 0.0) + g
                stats["kept_neighbor"] += 1
                continue

            steps = abs(dx) + abs(dy)
            if steps <= 0:
                continue
            r_step = e.r_ohm / float(steps)
            g_step = 1.0 / r_step
            cx, cy = x1, y1
            sx = 1 if dx > 0 else (-1 if dx < 0 else 0)
            sy = 1 if dy > 0 else (-1 if dy < 0 else 0)
            for _ in range(abs(dx)):
                nxp, nyp = cx + sx, cy
                a = _node_id(cx, cy, nx)
                b = _node_id(nxp, nyp, nx)
                if a > b:
                    a, b = b, a
                links_x[(a, b)] = links_x.get((a, b), 0.0) + g_step
                cx = nxp
            for _ in range(abs(dy)):
                nxp, nyp = cx, cy + sy
                a = _node_id(cx, cy, nx)
                b = _node_id(nxp, nyp, nx)
                if a > b:
                    a, b = b, a
                links_y[(a, b)] = links_y.get((a, b), 0.0) + g_step
                cy = nyp
            stats["split_non_neighbor"] += 1
        return links_x, links_y, stats

    m5_x, m5_y, stats_m5 = accumulate_links(edges_m5)
    m6_x, m6_y, stats_m6 = accumulate_links(edges_m6)

    merged_x = dict(m5_x)
    for k, g in m6_x.items():
        merged_x[k] = merged_x.get(k, 0.0) + g
    merged_y = dict(m5_y)
    for k, g in m6_y.items():
        merged_y[k] = merged_y.get(k, 0.0) + g

    rx = _estimate_pixel_axis_r(nx=nx, ny=ny, links=merged_x, axis="x")
    ry = _estimate_pixel_axis_r(nx=nx, ny=ny, links=merged_y, axis="y")

    via_count = [0] * (nx * ny)
    for _lineno, line in iter_netlist_lines(netlist):
        tokens = line.split()
        if not tokens or tokens[0][0] not in {"V", "v"}:
            continue
        if len(tokens) < 4:
            continue
        a, b = tokens[1], tokens[2]
        v = _safe_float(tokens[3])
        if v is None or abs(v) > 1e-15:
            continue
        if not (a.startswith("n") and b.startswith("n")):
            continue
        net_a = _net_id_from_node(a)
        net_b = _net_id_from_node(b)
        if net_a is None or net_b is None:
            continue
        if {net_a, net_b} != {1, 3}:
            continue
        xy_a = _parse_xy(a)
        xy_b = _parse_xy(b)
        if xy_a is None or xy_b is None:
            continue
        xa, ya = xy_a
        xb, yb = xy_b
        if xa != xb or ya != yb:
            continue
        if bbox is not None:
            bxmin, bxmax, bymin, bymax = bbox
            if not (bxmin <= xa <= bxmax and bymin <= ya <= bymax):
                continue
        px, py = _pixel_index(xa, ya, xmin=xmin, ymin=ymin, pitch=pitch)
        if not (0 <= px < nx and 0 <= py < ny):
            continue
        via_count[_node_id(px, py, nx)] += 1

    if rz_per_via <= 0:
        raise SystemExit("--rz-per-via must be > 0")
    rz = [float("inf")] * (nx * ny)
    for idx, c in enumerate(via_count):
        if c > 0:
            rz[idx] = rz_per_via / float(c)

    out_dir.mkdir(parents=True, exist_ok=True)
    meta = {
        "net": "VDD",
        "mode": "merged",
        "pitch": pitch,
        "origin": {"xmin": xmin, "ymin": ymin},
        "bbox": {"xmin": xmin, "xmax": xmax, "ymin": ymin, "ymax": ymax},
        "grid": {"nx": nx, "ny": ny},
        "rz_per_via": rz_per_via,
        "layers": {"m5": label_m5, "m6": label_m6},
        "edge_stats": {"m5": stats_m5, "m6": stats_m6},
        "notes": {
            "rx_ry": "Solved from merged neighbor link resistances (M5||M6).",
            "rz": "Computed from via_count between net 1 and net 3: rz = rz_per_via / via_count.",
        },
    }
    (out_dir / "meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    def write_matrix(path: Path, data: list[float]) -> None:
        with path.open("w", encoding="utf-8", newline="") as handle:
            w = csv.writer(handle)
            for py in range(ny):
                row = data[py * nx : (py + 1) * nx]
                w.writerow(row)

    def write_int_matrix(path: Path, data: list[int]) -> None:
        with path.open("w", encoding="utf-8", newline="") as handle:
            w = csv.writer(handle)
            for py in range(ny):
                row = data[py * nx : (py + 1) * nx]
                w.writerow(row)

    write_matrix(out_dir / "rx.csv", rx)
    write_matrix(out_dir / "ry.csv", ry)
    write_matrix(out_dir / "rz.csv", rz)
    write_int_matrix(out_dir / "via_count.csv", via_count)


def plot_pixel_csv_heatmap(in_csv: Path, out_png: Path, *, title: str) -> None:
    _ensure_mpl_cache_dir()
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt  # noqa: E402

    rows: list[list[float]] = []
    with in_csv.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            parts = [p.strip() for p in line.split(",")]
            row: list[float] = []
            for p in parts:
                try:
                    row.append(float(p))
                except ValueError:
                    row.append(float("nan"))
            rows.append(row)
    if not rows:
        raise SystemExit(f"Empty csv: {in_csv}")

    fig, ax = plt.subplots(figsize=(7, 6), dpi=200)
    im = ax.imshow(rows, origin="lower", aspect="auto")
    ax.set_title(title)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)

def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Parse and visualize per-layer resistor networks from a SPICE PDN netlist.")
    parser.add_argument("--netlist", type=Path, help="Path to SPICE netlist (e.g. ibmpg1.spice)")
    parser.add_argument("--plot-branches", type=Path, help="Plot a PDNNetwork branches CSV (layer,row,col,direction,resistance,node1,node2)")
    parser.add_argument("--out-branches-png", type=Path, help="Output PNG path for --plot-branches")
    parser.add_argument("--show-branches", action="store_true", help="Show an interactive Matplotlib window for --plot-branches")
    parser.add_argument("--branches-stride", type=int, default=1, help="Plot every Nth branch (reduces clutter)")
    parser.add_argument("--branches-lw", type=float, default=0.6, help="Line width for branches plot (Matplotlib)")
    parser.add_argument("--branches-alpha", type=float, default=0.35, help="Line alpha for branches plot (Matplotlib)")
    parser.add_argument("--branches-dpi", type=int, default=220, help="PNG DPI for branches plot (Matplotlib)")
    parser.add_argument("--list-layers", action="store_true", help="List detected '* layer:' sections and resistor counts")
    parser.add_argument("--layer", help="Layer selector (substring match against the '* layer:' label)")
    parser.add_argument("--all", action="store_true", help="Render all detected layers")
    parser.add_argument("--out", type=Path, help="Output SVG path (single layer mode)")
    parser.add_argument("--outdir", type=Path, help="Output directory (all-layers mode)")
    parser.add_argument("--out-3d-png", type=Path, help="Output PNG path (Matplotlib 3D line plot; includes all layers unless --layer is set)")
    parser.add_argument("--z-step", type=float, default=2000.0, help="Z spacing between metal layers (3D export)")
    parser.add_argument("--no-vias", action="store_true", help="Do not include via shorts in 3D export")
    parser.add_argument("--show-3d", action="store_true", help="Show an interactive 3D Matplotlib window (rotate with mouse)")
    parser.add_argument("--3d-elev", dest="d3_elev", type=float, default=25.0, help="3D view elevation (Matplotlib)")
    parser.add_argument("--3d-azim", dest="d3_azim", type=float, default=-60.0, help="3D view azimuth (Matplotlib)")
    parser.add_argument("--3d-stride", dest="d3_stride", type=int, default=1, help="Plot every Nth metal segment in 3D (speed/size tradeoff)")
    parser.add_argument("--3d-lw", dest="d3_lw", type=float, default=0.6, help="3D line width (Matplotlib)")
    parser.add_argument("--3d-alpha", dest="d3_alpha", type=float, default=0.65, help="3D line alpha (Matplotlib)")
    parser.add_argument("--3d-dpi", dest="d3_dpi", type=int, default=220, help="3D PNG DPI (Matplotlib)")
    parser.add_argument("--vias-stride", type=int, default=25, help="Plot every Nth via in 3D (reduces clutter)")
    parser.add_argument("--vias-lw", type=float, default=0.25, help="Via line width in 3D (Matplotlib)")
    parser.add_argument("--vias-alpha", type=float, default=0.08, help="Via alpha in 3D (Matplotlib)")
    parser.add_argument("--3d-stats", dest="d3_stats", action="store_true", help="Print how many metal/via segments are plotted in 3D")
    parser.add_argument("--export-csv", type=Path, help="Optional directory to write per-layer edge CSV files")
    parser.add_argument("--bbox", nargs=4, metavar=("XMIN", "XMAX", "YMIN", "YMAX"), help="Crop to a coordinate window")
    parser.add_argument("--min-r", type=float, default=None, help="Filter: minimum resistor value (ohm)")
    parser.add_argument("--max-r", type=float, default=None, help="Filter: maximum resistor value (ohm)")
    parser.add_argument("--max-edges", type=int, default=None, help="Stop after collecting this many edges (for huge nets)")
    parser.add_argument("--width", type=int, default=1400, help="SVG width in pixels")
    parser.add_argument("--height", type=int, default=1400, help="SVG height in pixels")
    parser.add_argument("--margin", type=int, default=20, help="SVG margin in pixels")
    parser.add_argument("--export-pixel-vdd", type=Path, help="Export a VDD pixel model (Rx/Ry per layer + Rz via) to this directory")
    parser.add_argument("--export-pixel-vdd-merged", type=Path, help="Export a merged single-layer VDD pixel model (rx/ry/rz) to this directory")
    parser.add_argument("--pitch", type=float, default=None, help="Pixel pitch (same units as node coordinates); default auto from bbox")
    parser.add_argument("--target-pixels", type=int, default=64, help="Auto-pitch targets this many pixels across max dimension")
    parser.add_argument("--rz-per-via", type=float, default=1e-4, help="Assumed resistance per single via (ohm) when computing Rz = rz_per_via / via_count")
    parser.add_argument("--plot-pixel", type=Path, help="Given an export directory, write rx/ry/rz heatmaps (PNG) into it")
    args = parser.parse_args(argv)

    if args.plot_branches is not None:
        show = args.show_branches or args.out_branches_png is None
        plot_branches_virtual_grid(
            branches_csv=args.plot_branches,
            out_png=args.out_branches_png,
            show=show,
            stride=args.branches_stride,
            lw=args.branches_lw,
            alpha=args.branches_alpha,
            dpi=args.branches_dpi,
        )
        if args.out_branches_png is not None:
            print(f"Wrote {args.out_branches_png}")
        return 0

    needs_netlist = any(
        [
            args.list_layers,
            args.all,
            args.layer is not None,
            args.out is not None,
            args.outdir is not None,
            args.out_3d_png is not None,
            args.show_3d,
            args.export_csv is not None,
            args.export_pixel_vdd is not None,
            args.export_pixel_vdd_merged is not None,
        ]
    )
    if needs_netlist:
        if args.netlist is None:
            raise SystemExit("--netlist is required for this operation")
        if not args.netlist.exists():
            raise SystemExit(f"Netlist not found: {args.netlist}")

    bbox = _parse_bbox(args.bbox) if args.bbox is not None else None

    if args.export_pixel_vdd is not None:
        export_pixel_model_vdd(
            args.netlist,
            args.export_pixel_vdd,
            pitch=args.pitch,
            target_pixels=args.target_pixels,
            rz_per_via=args.rz_per_via,
            bbox=bbox,
            max_edges=args.max_edges,
        )
        print(f"Wrote VDD pixel model to {args.export_pixel_vdd}")
        return 0

    if args.export_pixel_vdd_merged is not None:
        export_pixel_model_vdd_merged(
            args.netlist,
            args.export_pixel_vdd_merged,
            pitch=args.pitch,
            target_pixels=args.target_pixels,
            rz_per_via=args.rz_per_via,
            bbox=bbox,
            max_edges=args.max_edges,
        )
        print(f"Wrote merged VDD pixel model to {args.export_pixel_vdd_merged}")
        return 0

    if args.plot_pixel is not None:
        # Heuristic: if merged export exists, plot rx/ry/rz. Otherwise plot per-layer exports if present.
        rx = args.plot_pixel / "rx.csv"
        ry = args.plot_pixel / "ry.csv"
        rz = args.plot_pixel / "rz.csv"
        if rx.exists() and ry.exists() and rz.exists():
            plot_pixel_csv_heatmap(rx, args.plot_pixel / "rx.png", title="VDD merged Rx")
            plot_pixel_csv_heatmap(ry, args.plot_pixel / "ry.png", title="VDD merged Ry")
            plot_pixel_csv_heatmap(rz, args.plot_pixel / "rz.png", title="VDD merged Rz (via)")
            print(f"Wrote heatmaps to {args.plot_pixel}")
            return 0
        raise SystemExit(f"Did not find merged rx/ry/rz csv in {args.plot_pixel}")

    did_any_3d = False
    if args.out_3d_png is not None or args.show_3d:
        segments = extract_3d_segments(
            args.netlist,
            layer_selector=args.layer,
            bbox=bbox,
            min_r=args.min_r,
            max_r=args.max_r,
            max_edges=args.max_edges,
            include_vias=not args.no_vias,
            z_step=args.z_step,
        )
        if args.d3_stats:
            stats = segment_stats(segments)
            print(f"3D segments: total={stats['total']} metal={stats['metal']} vias={stats['vias']}")
        if args.out_3d_png is not None:
            write_3d_png(
                args.out_3d_png,
                segments=segments,
                elev=args.d3_elev,
                azim=args.d3_azim,
                stride=args.d3_stride,
                linewidth=args.d3_lw,
                alpha=args.d3_alpha,
                vias_stride=args.vias_stride,
                vias_linewidth=args.vias_lw,
                vias_alpha=args.vias_alpha,
                dpi=args.d3_dpi,
                figsize=(10.0, 10.0),
            )
            print(f"Wrote {args.out_3d_png}")
            did_any_3d = True
        if args.show_3d:
            fig, _ax = render_3d_matplotlib(
                segments=segments,
                elev=args.d3_elev,
                azim=args.d3_azim,
                stride=args.d3_stride,
                linewidth=args.d3_lw,
                alpha=args.d3_alpha,
                vias_stride=args.vias_stride,
                vias_linewidth=args.vias_lw,
                vias_alpha=args.vias_alpha,
                figsize=(10.0, 10.0),
                dpi=args.d3_dpi,
                interactive=True,
            )
            from matplotlib import pyplot as plt

            did_any_3d = True
            plt.show()
            plt.close(fig)

        if did_any_3d and not args.all and args.out is None and args.outdir is None:
            return 0

    if args.list_layers:
        layers = scan_layers(args.netlist)
        if not layers:
            print("No '* layer:' sections found.")
            return 0
        print("Detected layers:")
        for s in layers:
            print(f"- line {s.start_line:>6}: {s.label}  (resistors: {s.resistor_edges})")
        return 0

    if args.all:
        layers = scan_layers(args.netlist)
        if not layers:
            raise SystemExit("No '* layer:' sections found; cannot --all render.")
        if args.outdir is None:
            raise SystemExit("--all requires --outdir")
        for s in layers:
            selected_label, edges = extract_layer_edges(
                args.netlist,
                layer_wanted=s.label,
                bbox=bbox,
                min_r=args.min_r,
                max_r=args.max_r,
                max_edges=args.max_edges,
            )
            out_path = args.outdir / f"{_slug(selected_label)}.svg"
            write_svg(
                out_path,
                title=f"{selected_label} | edges={len(edges)}",
                edges=edges,
                width_px=args.width,
                height_px=args.height,
                margin_px=args.margin,
            )
            if args.export_csv is not None:
                export_csv(args.export_csv, label=selected_label, edges=edges)
            print(f"Wrote {out_path}")
        return 0

    if not args.layer:
        raise SystemExit("Provide --layer (or use --list-layers / --all).")
    if args.out is None:
        raise SystemExit("Single layer mode requires --out.")

    selected_label, edges = extract_layer_edges(
        args.netlist,
        layer_wanted=args.layer,
        bbox=bbox,
        min_r=args.min_r,
        max_r=args.max_r,
        max_edges=args.max_edges,
    )
    write_svg(
        args.out,
        title=f"{selected_label} | edges={len(edges)}",
        edges=edges,
        width_px=args.width,
        height_px=args.height,
        margin_px=args.margin,
    )
    if args.export_csv is not None:
        export_csv(args.export_csv, label=selected_label, edges=edges)
    print(f"Wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
