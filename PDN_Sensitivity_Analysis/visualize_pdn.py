#!/usr/bin/env python3
"""
PDN Network Visualization Tool

This script visualizes the PDN network structure, showing:
- Grid nodes and their connections
- Resistance values on branches
- Network topology
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import sys
import os

class PDNVisualizer:
    def __init__(self, rows, cols, layers=1):
        self.rows = rows
        self.cols = cols
        self.layers = layers
        
    def visualize_grid_2d(self, pixel_models=None, save_file=None, show=True):
        """
        Visualize 2D PDN grid (top view)
        
        Args:
            pixel_models: Dictionary of (row, col) -> (Rx, Ry, Rz) or None
            save_file: Filename to save the plot (optional)
            show: Whether to display the plot
        """
        fig, ax = plt.subplots(1, 1, figsize=(12, 10))
        
        # Draw grid
        node_positions = {}
        node_id = 0
        
        # Calculate node positions
        for r in range(self.rows):
            for c in range(self.cols):
                x = c
                y = self.rows - 1 - r  # Flip Y axis for better visualization
                node_positions[(r, c)] = (x, y)
                
                # Draw node
                circle = plt.Circle((x, y), 0.1, color='blue', zorder=3)
                ax.add_patch(circle)
                ax.text(x, y, f'({r},{c})', fontsize=6, ha='center', va='center', 
                       color='white', weight='bold', zorder=4)
        
        # Draw horizontal branches (x-direction)
        for r in range(self.rows):
            for c in range(self.cols - 1):
                pos1 = node_positions[(r, c)]
                pos2 = node_positions[(r, c + 1)]
                
                # Get resistance value
                if pixel_models and (r, c) in pixel_models:
                    rx = pixel_models[(r, c)][0]
                    color = self._get_resistance_color(rx, 0.05, 0.2)
                    label = f'Rx={rx:.3f}'
                else:
                    color = 'gray'
                    label = 'Rx'
                
                # Draw branch
                ax.plot([pos1[0], pos2[0]], [pos1[1], pos2[1]], 
                       color=color, linewidth=2, zorder=1)
                
                # Add label
                mid_x = (pos1[0] + pos2[0]) / 2
                mid_y = (pos1[1] + pos2[1]) / 2
                ax.text(mid_x, mid_y + 0.15, label, fontsize=7, 
                       ha='center', va='bottom', color=color, weight='bold')
        
        # Draw vertical branches (y-direction)
        for r in range(self.rows - 1):
            for c in range(self.cols):
                pos1 = node_positions[(r, c)]
                pos2 = node_positions[(r + 1, c)]
                
                # Get resistance value
                if pixel_models and (r, c) in pixel_models:
                    ry = pixel_models[(r, c)][1]
                    color = self._get_resistance_color(ry, 0.05, 0.2)
                    label = f'Ry={ry:.3f}'
                else:
                    color = 'gray'
                    label = 'Ry'
                
                # Draw branch
                ax.plot([pos1[0], pos2[0]], [pos1[1], pos2[1]], 
                       color=color, linewidth=2, zorder=1)
                
                # Add label
                mid_x = (pos1[0] + pos2[0]) / 2
                mid_y = (pos1[1] + pos2[1]) / 2
                ax.text(mid_x + 0.15, mid_y, label, fontsize=7, 
                       ha='left', va='center', color=color, weight='bold')
        
        # Set up axes
        ax.set_xlim(-0.5, self.cols - 0.5)
        ax.set_ylim(-0.5, self.rows - 0.5)
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.set_xlabel('Column', fontsize=12)
        ax.set_ylabel('Row', fontsize=12)
        ax.set_title(f'PDN Network Grid ({self.rows}×{self.cols})', fontsize=14, weight='bold')
        
        # Add legend
        legend_elements = [
            mpatches.Patch(color='blue', label='Grid Node'),
            mpatches.Patch(color='green', label='Low Resistance'),
            mpatches.Patch(color='orange', label='Medium Resistance'),
            mpatches.Patch(color='red', label='High Resistance'),
        ]
        ax.legend(handles=legend_elements, loc='upper right')
        
        plt.tight_layout()
        
        if save_file:
            plt.savefig(save_file, dpi=300, bbox_inches='tight')
            print(f"Visualization saved to {save_file}")
        
        if show:
            plt.show()
        else:
            plt.close()
    
    def visualize_resistance_heatmap(self, pixel_models, direction='x', save_file=None, show=True):
        """
        Visualize resistance values as a heatmap
        
        Args:
            pixel_models: Dictionary of (row, col) -> (Rx, Ry, Rz)
            direction: 'x', 'y', or 'z' to show which resistance
            save_file: Filename to save the plot
            show: Whether to display the plot
        """
        fig, ax = plt.subplots(1, 1, figsize=(10, 8))
        
        # Create resistance matrix
        resistance_matrix = np.zeros((self.rows, self.cols))
        
        for r in range(self.rows):
            for c in range(self.cols):
                if (r, c) in pixel_models:
                    if direction == 'x':
                        resistance_matrix[r, c] = pixel_models[(r, c)][0]
                    elif direction == 'y':
                        resistance_matrix[r, c] = pixel_models[(r, c)][1]
                    else:
                        resistance_matrix[r, c] = pixel_models[(r, c)][2]
        
        # Create heatmap
        im = ax.imshow(resistance_matrix, cmap='YlOrRd', aspect='auto', origin='upper')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label(f'R{direction} Resistance (Ω)', fontsize=12)
        
        # Add text annotations
        for r in range(self.rows):
            for c in range(self.cols):
                text = ax.text(c, r, f'{resistance_matrix[r, c]:.3f}',
                             ha="center", va="center", color="black", fontsize=8)
        
        ax.set_xlabel('Column', fontsize=12)
        ax.set_ylabel('Row', fontsize=12)
        ax.set_title(f'R{direction.upper()} Resistance Heatmap', fontsize=14, weight='bold')
        ax.set_xticks(range(self.cols))
        ax.set_yticks(range(self.rows))
        
        plt.tight_layout()
        
        if save_file:
            plt.savefig(save_file, dpi=300, bbox_inches='tight')
            print(f"Heatmap saved to {save_file}")
        
        if show:
            plt.show()
        else:
            plt.close()
    
    def visualize_network_graph(self, branches=None, save_file=None, show=True):
        """
        Visualize network as a graph (using networkx-style visualization)
        
        Args:
            branches: List of (node1, node2, resistance) tuples
            save_file: Filename to save the plot
            show: Whether to display the plot
        """
        fig, ax = plt.subplots(1, 1, figsize=(12, 10))
        
        # Calculate node positions
        node_positions = {}
        for r in range(self.rows):
            for c in range(self.cols):
                x = c * 2
                y = (self.rows - 1 - r) * 2
                node_positions[(r, c)] = (x, y)
        
        # Draw branches
        if branches:
            for node1, node2, resistance in branches:
                # Convert node IDs to positions (simplified)
                # In real implementation, you'd map node IDs to (r, c)
                pass
        
        # Draw all possible branches
        for r in range(self.rows):
            for c in range(self.cols):
                pos = node_positions[(r, c)]
                
                # Draw node
                circle = plt.Circle(pos, 0.15, color='blue', zorder=3)
                ax.add_patch(circle)
                ax.text(pos[0], pos[1], f'{r},{c}', fontsize=7, 
                       ha='center', va='center', color='white', weight='bold', zorder=4)
                
                # Draw horizontal connection
                if c < self.cols - 1:
                    pos2 = node_positions[(r, c + 1)]
                    ax.plot([pos[0], pos2[0]], [pos[1], pos2[1]], 
                           'k-', linewidth=1, alpha=0.5, zorder=1)
                
                # Draw vertical connection
                if r < self.rows - 1:
                    pos2 = node_positions[(r + 1, c)]
                    ax.plot([pos[0], pos2[0]], [pos[1], pos2[1]], 
                           'k-', linewidth=1, alpha=0.5, zorder=1)
        
        ax.set_aspect('equal')
        ax.set_xlim(-0.5, (self.cols - 1) * 2 + 0.5)
        ax.set_ylim(-0.5, (self.rows - 1) * 2 + 0.5)
        ax.set_title(f'PDN Network Graph ({self.rows}×{self.cols})', fontsize=14, weight='bold')
        ax.axis('off')
        
        plt.tight_layout()
        
        if save_file:
            plt.savefig(save_file, dpi=300, bbox_inches='tight')
            print(f"Graph saved to {save_file}")
        
        if show:
            plt.show()
        else:
            plt.close()
    
    def _get_resistance_color(self, resistance, min_r, max_r):
        """Get color based on resistance value"""
        normalized = (resistance - min_r) / (max_r - min_r) if max_r > min_r else 0.5
        normalized = max(0, min(1, normalized))  # Clamp to [0, 1]
        
        if normalized < 0.33:
            return 'green'  # Low resistance
        elif normalized < 0.67:
            return 'orange'  # Medium resistance
        else:
            return 'red'  # High resistance


def create_sample_pdn(rows=5, cols=5):
    """Create a sample PDN with varying resistances"""
    pixel_models = {}
    base_rx, base_ry, base_rz = 0.1, 0.1, 0.5
    
    for r in range(rows):
        for c in range(cols):
            # Vary resistance slightly for visualization
            rx = base_rx * (1 + 0.1 * (r + c) / (rows + cols))
            ry = base_ry * (1 + 0.1 * (r + c) / (rows + cols))
            rz = base_rz
            pixel_models[(r, c)] = (rx, ry, rz)
    
    return pixel_models


def main():
    """Main function for command-line usage"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Visualize PDN Network')
    parser.add_argument('--rows', type=int, default=5, help='Number of rows')
    parser.add_argument('--cols', type=int, default=5, help='Number of columns')
    parser.add_argument('--layers', type=int, default=1, help='Number of layers')
    parser.add_argument('--type', choices=['grid', 'heatmap', 'graph'], 
                       default='grid', help='Visualization type')
    parser.add_argument('--direction', choices=['x', 'y', 'z'], 
                       default='x', help='Resistance direction for heatmap')
    parser.add_argument('--output', type=str, help='Output filename (optional)')
    parser.add_argument('--no-show', action='store_true', 
                       help='Do not display plot (useful when saving)')
    
    args = parser.parse_args()
    
    # Create visualizer
    viz = PDNVisualizer(args.rows, args.cols, args.layers)
    
    # Create sample data
    pixel_models = create_sample_pdn(args.rows, args.cols)
    
    # Generate visualization
    if args.type == 'grid':
        viz.visualize_grid_2d(pixel_models, save_file=args.output, 
                            show=not args.no_show)
    elif args.type == 'heatmap':
        viz.visualize_resistance_heatmap(pixel_models, direction=args.direction,
                                       save_file=args.output, show=not args.no_show)
    elif args.type == 'graph':
        viz.visualize_network_graph(save_file=args.output, show=not args.no_show)


if __name__ == '__main__':
    main()

