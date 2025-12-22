#!/usr/bin/env python3
"""
Visualize PDN network with voltage sources and current loads highlighted.

Usage:
    python3 tools/visualize_pdn_loads.py --branches out_voltspot/pdn_with_loads_branches.csv --nodes out_voltspot/pdn_with_loads_branches_nodes.csv --out out_voltspot/pdn_visualization.png
"""

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, Set, Tuple

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


def read_branches_csv(branches_path: Path) -> list:
    """Read branches CSV file."""
    branches = []
    with open(branches_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            branches.append({
                'layer': int(row['layer']),
                'row': int(row['row']),
                'col': int(row['col']),
                'direction': row['direction'],
                'resistance': float(row['resistance']),
                'node1': int(row['node1']),
                'node2': int(row['node2']),
            })
    return branches


def read_nodes_csv(nodes_path: Path) -> Dict[Tuple[int, int, int], dict]:
    """Read nodes CSV file and return a map from (layer, row, col) to node info."""
    nodes = {}
    with open(nodes_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = (int(row['layer']), int(row['row']), int(row['col']))
            nodes[key] = {
                'pgNet': row['pgNet'],
                'nodeId': int(row['nodeId']),
                'voltage': float(row['voltage']),
                'currentLoad': float(row['currentLoad']),
                'isPad': int(row.get('isPad', '0')),
            }
    return nodes


def visualize_pdn(branches: list, nodes: Dict[Tuple[int, int, int], dict],
                  out_path: Path, layer: int = 0, stride: int = 1,
                  show_branches: bool = True, show_nodes: bool = True):
    """Visualize PDN network with voltage sources and current loads."""
    
    # Filter branches for the specified layer
    layer_branches = [b for b in branches if b['layer'] == layer]
    
    # Get grid dimensions
    max_row = max(b['row'] for b in layer_branches) if layer_branches else 0
    max_col = max(b['col'] for b in layer_branches) if layer_branches else 0
    
    # Create figure
    fig, ax = plt.subplots(figsize=(16, 12))
    
    # Draw branches (as lines)
    if show_branches:
        for branch in layer_branches[::stride]:
            row = branch['row']
            col = branch['col']
            direction = branch['direction']
            
            if direction == 'x':
                # Horizontal branch
                ax.plot([col, col + 1], [row, row], 'b-', alpha=0.3, linewidth=0.5)
            elif direction == 'y':
                # Vertical branch
                ax.plot([col, col], [row, row + 1], 'b-', alpha=0.3, linewidth=0.5)
    
    # Collect node positions and types
    voltage_source_nodes = []
    current_load_nodes = []
    regular_nodes = []
    
    for (l, r, c), node_info in nodes.items():
        if l != layer:
            continue
        
        # Prefer explicit pad flag (so GND pads at 0V are still treated as voltage sources)
        if node_info.get('isPad', 0) != 0:
            voltage_source_nodes.append((c, r, node_info['voltage'], node_info['currentLoad']))
        elif node_info['currentLoad'] != 0.0:
            current_load_nodes.append((c, r, node_info['voltage'], node_info['currentLoad']))
        else:
            regular_nodes.append((c, r))
    
    # Draw nodes
    if show_nodes:
        # Regular nodes (small gray dots)
        if regular_nodes:
            x, y = zip(*regular_nodes)
            ax.scatter(x, y, c='gray', s=2, alpha=0.3, marker='o', label='Regular nodes')
        
        # Current load nodes (colored by load magnitude)
        if current_load_nodes:
            x, y, v, loads = zip(*current_load_nodes)
            loads_array = np.array(loads)
            # Normalize for color mapping
            if loads_array.max() > 0:
                normalized_loads = loads_array / loads_array.max()
            else:
                normalized_loads = loads_array
            
            scatter = ax.scatter(x, y, c=normalized_loads, s=50, alpha=0.7, 
                               marker='s', cmap='YlOrRd', vmin=0, vmax=1,
                               label='Current load nodes', edgecolors='black', linewidths=0.5)
            
            # Add colorbar for current loads
            cbar = plt.colorbar(scatter, ax=ax, label='Normalized Current Load')
        
        # Voltage source nodes (red stars)
        if voltage_source_nodes:
            x, y, v, loads = zip(*voltage_source_nodes)
            # Size based on voltage magnitude
            sizes = [abs(voltage) * 100 for voltage in v]
            ax.scatter(x, y, c='red', s=sizes, alpha=0.8, marker='*', 
                      edgecolors='darkred', linewidths=1.5, 
                      label='Voltage source nodes', zorder=10)
    
    # Set labels and title
    ax.set_xlabel('Column', fontsize=12)
    ax.set_ylabel('Row', fontsize=12)
    ax.set_title(f'PDN Network Visualization (Layer {layer})\n'
                f'Voltage Sources: {len(voltage_source_nodes)}, '
                f'Current Loads: {len(current_load_nodes)}', fontsize=14)
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal', adjustable='box')
    
    # NOTE: CSV rows are already in VoltSpot coordinate convention (bottom-origin),
    # so we do NOT invert y-axis here.
    
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"Visualization saved to {out_path}")
    
    # Print statistics
    print(f"\nStatistics:")
    print(f"  Total nodes: {len([n for k, n in nodes.items() if k[0] == layer])}")
    print(f"  Voltage source nodes: {len(voltage_source_nodes)}")
    print(f"  Current load nodes: {len(current_load_nodes)}")
    print(f"  Regular nodes: {len(regular_nodes)}")
    
    if voltage_source_nodes:
        voltages = [v for _, _, v, _ in voltage_source_nodes]
        print(f"  Voltage range: {min(voltages):.3f}V to {max(voltages):.3f}V")
    
    if current_load_nodes:
        loads = [l for _, _, _, l in current_load_nodes]
        print(f"  Current load range: {min(loads):.6f}A to {max(loads):.6f}A")
        print(f"  Total current load: {sum(loads):.6f}A")


def main():
    parser = argparse.ArgumentParser(
        description='Visualize PDN network with voltage sources and current loads'
    )
    parser.add_argument('--branches', type=Path, required=True,
                       help='Branches CSV file (e.g., pdn_with_loads_branches.csv)')
    parser.add_argument('--nodes', type=Path, required=True,
                       help='Nodes CSV file (e.g., pdn_with_loads_branches_nodes.csv)')
    parser.add_argument('--out', type=Path, required=True,
                       help='Output image file (PNG)')
    parser.add_argument('--layer', type=int, default=0,
                       help='Layer to visualize (default: 0)')
    parser.add_argument('--stride', type=int, default=1,
                       help='Plot every Nth branch (default: 1, use larger values for dense networks)')
    parser.add_argument('--no-branches', action='store_true',
                       help='Hide branch lines')
    parser.add_argument('--no-nodes', action='store_true',
                       help='Hide node markers')
    
    args = parser.parse_args()
    
    # Check files exist
    if not args.branches.exists():
        print(f"Error: Branches file not found: {args.branches}", file=sys.stderr)
        sys.exit(1)
    
    if not args.nodes.exists():
        print(f"Error: Nodes file not found: {args.nodes}", file=sys.stderr)
        sys.exit(1)
    
    # Read data
    print(f"Reading branches from {args.branches}...")
    branches = read_branches_csv(args.branches)
    print(f"  Found {len(branches)} branches")
    
    print(f"Reading nodes from {args.nodes}...")
    nodes = read_nodes_csv(args.nodes)
    print(f"  Found {len(nodes)} nodes")
    
    # Create output directory if needed
    args.out.parent.mkdir(parents=True, exist_ok=True)
    
    # Visualize
    print(f"\nGenerating visualization...")
    visualize_pdn(branches, nodes, args.out, layer=args.layer, stride=args.stride,
                  show_branches=not args.no_branches, show_nodes=not args.no_nodes)
    
    print("\nDone!")


if __name__ == '__main__':
    main()

