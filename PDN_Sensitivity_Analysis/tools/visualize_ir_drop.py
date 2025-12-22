#!/usr/bin/env python3
"""
Visualize IR drop from gridIR file (VoltSpot format).

Usage:
    python3 tools/visualize_ir_drop.py --gridir out_voltspot/ir_drop.gridIR --out out_voltspot/ir_drop_visualization.png
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np


def parse_gridir_file(gridir_path: Path) -> Dict[int, np.ndarray]:
    """
    Parse gridIR file and return a dictionary mapping layer -> 2D array of IR drop percentages.
    
    Format:
        #Layer:N
        col row drop_ratio%
        ...
    """
    ir_drop_data = {}
    current_layer = None
    current_data = []
    max_row = 0
    max_col = 0
    
    with open(gridir_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                # Empty line: save current layer data
                if current_layer is not None and current_data:
                    # Find dimensions
                    rows = max(r for _, r, _ in current_data) + 1
                    cols = max(c for c, _, _ in current_data) + 1
                    
                    # Create 2D array
                    layer_array = np.zeros((rows, cols))
                    for col, row, drop in current_data:
                        layer_array[row, col] = drop
                    
                    ir_drop_data[current_layer] = layer_array
                    current_data = []
                    current_layer = None
                continue
            
            if line.startswith('#Layer:'):
                # New layer
                current_layer = int(line.split(':')[1])
                current_data = []
            else:
                # Data line: col row drop_ratio%
                parts = line.split()
                if len(parts) >= 3:
                    col = int(parts[0])
                    row = int(parts[1])
                    drop = float(parts[2])
                    current_data.append((col, row, drop))
                    max_row = max(max_row, row)
                    max_col = max(max_col, col)
    
    # Handle last layer if file doesn't end with empty line
    if current_layer is not None and current_data:
        rows = max(r for _, r, _ in current_data) + 1
        cols = max(c for c, _, _ in current_data) + 1
        
        layer_array = np.zeros((rows, cols))
        for col, row, drop in current_data:
            layer_array[row, col] = drop
        
        ir_drop_data[current_layer] = layer_array
    
    return ir_drop_data


def visualize_ir_drop(ir_drop_data: Dict[int, np.ndarray], out_path: Path,
                      layer: int = None, cmap: str = 'hot', 
                      vmin: float = None, vmax: float = None,
                      show_colorbar: bool = True):
    """
    Visualize IR drop data as heatmaps.
    
    Args:
        ir_drop_data: Dictionary mapping layer -> 2D array of IR drop percentages
        out_path: Output image path
        layer: Specific layer to visualize (None = all layers in subplots)
        cmap: Colormap name
        vmin: Minimum value for colormap (None = auto)
        vmax: Maximum value for colormap (None = auto)
        show_colorbar: Whether to show colorbar
    """
    if not ir_drop_data:
        print("Error: No IR drop data to visualize", file=sys.stderr)
        return
    
    # Determine value range
    all_values = []
    for layer_data in ir_drop_data.values():
        all_values.extend(layer_data.flatten())
    
    if vmin is None:
        vmin = min(all_values) if all_values else 0.0
    if vmax is None:
        vmax = max(all_values) if all_values else 1.0
    
    # Determine layout
    if layer is not None:
        # Single layer
        if layer not in ir_drop_data:
            print(f"Error: Layer {layer} not found in data", file=sys.stderr)
            return
        
        fig, ax = plt.subplots(figsize=(12, 10))
        data = ir_drop_data[layer]
        
        im = ax.imshow(data, cmap=cmap, vmin=vmin, vmax=vmax, origin='lower', aspect='auto')
        ax.set_title(f'IR Drop - Layer {layer}\nMax: {np.max(data):.3f}%, Avg: {np.mean(data):.3f}%', 
                    fontsize=14, fontweight='bold')
        ax.set_xlabel('Column', fontsize=12)
        ax.set_ylabel('Row', fontsize=12)
        
        if show_colorbar:
            cbar = plt.colorbar(im, ax=ax)
            cbar.set_label('IR Drop (%)', fontsize=12)
        
    else:
        # All layers in subplots
        num_layers = len(ir_drop_data)
        layers = sorted(ir_drop_data.keys())
        
        # Calculate grid layout
        cols = min(3, num_layers)
        rows = (num_layers + cols - 1) // cols
        
        fig, axes = plt.subplots(rows, cols, figsize=(6*cols, 5*rows))
        if num_layers == 1:
            axes = [axes]
        else:
            axes = axes.flatten()
        
        for idx, layer_num in enumerate(layers):
            ax = axes[idx]
            data = ir_drop_data[layer_num]
            
            im = ax.imshow(data, cmap=cmap, vmin=vmin, vmax=vmax, origin='lower', aspect='auto')
            ax.set_title(f'Layer {layer_num}\nMax: {np.max(data):.3f}%, Avg: {np.mean(data):.3f}%',
                        fontsize=11)
            ax.set_xlabel('Column', fontsize=10)
            ax.set_ylabel('Row', fontsize=10)
            
            if show_colorbar and idx == 0:
                # Add colorbar to first subplot
                cbar = plt.colorbar(im, ax=ax)
                cbar.set_label('IR Drop (%)', fontsize=10)
        
        # Hide unused subplots
        for idx in range(num_layers, len(axes)):
            axes[idx].axis('off')
        
        fig.suptitle('IR Drop Visualization', fontsize=16, fontweight='bold', y=0.995)
    
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"IR drop visualization saved to: {out_path}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Visualize IR drop from gridIR file (VoltSpot format)'
    )
    parser.add_argument('--gridir', type=Path, required=True,
                       help='Input gridIR file (e.g., ir_drop.gridIR)')
    parser.add_argument('--out', type=Path, required=True,
                       help='Output image file (PNG)')
    parser.add_argument('--layer', type=int, default=None,
                       help='Specific layer to visualize (default: all layers)')
    parser.add_argument('--cmap', type=str, default='hot',
                       help='Colormap name (default: hot, options: hot, cool, viridis, plasma, etc.)')
    parser.add_argument('--vmin', type=float, default=None,
                       help='Minimum IR drop value for colormap (default: auto)')
    parser.add_argument('--vmax', type=float, default=None,
                       help='Maximum IR drop value for colormap (default: auto)')
    parser.add_argument('--no-colorbar', action='store_true',
                       help='Hide colorbar')
    
    args = parser.parse_args()
    
    # Check file exists
    if not args.gridir.exists():
        print(f"Error: gridIR file not found: {args.gridir}", file=sys.stderr)
        sys.exit(1)
    
    # Parse gridIR file
    print(f"Reading gridIR file: {args.gridir}...")
    ir_drop_data = parse_gridir_file(args.gridir)
    
    if not ir_drop_data:
        print("Error: No data found in gridIR file", file=sys.stderr)
        sys.exit(1)
    
    print(f"  Found {len(ir_drop_data)} layer(s): {sorted(ir_drop_data.keys())}")
    for layer_num, data in ir_drop_data.items():
        print(f"    Layer {layer_num}: {data.shape[0]}x{data.shape[1]} grid, "
              f"IR drop range: {np.min(data):.3f}% - {np.max(data):.3f}%, "
              f"avg: {np.mean(data):.3f}%")
    
    # Visualize
    print(f"\nGenerating visualization...")
    visualize_ir_drop(ir_drop_data, args.out, layer=args.layer, 
                     cmap=args.cmap, vmin=args.vmin, vmax=args.vmax,
                     show_colorbar=not args.no_colorbar)
    
    print("Done!")


if __name__ == '__main__':
    main()

