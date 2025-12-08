#!/usr/bin/env python3
"""
Plot and Analyze Sensitivity Analysis Results (No pandas dependency)

This script reads CSV files from sensitivity analysis and creates visualizations
using only standard library and matplotlib.
"""

import csv
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

# Set style for better-looking plots
try:
    plt.style.use('seaborn-v0_8-darkgrid')
except:
    plt.style.use('default')
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 10

class SensitivityPlotter:
    def __init__(self, output_dir='.'):
        self.output_dir = output_dir
        self.results = {}
        
    def load_csv(self, filename):
        """Load CSV file and return data as dictionary"""
        if not os.path.exists(filename):
            print(f"Warning: {filename} not found, skipping...")
            return None
        
        data = {'headers': [], 'rows': []}
        
        with open(filename, 'r') as f:
            reader = csv.reader(f)
            data['headers'] = next(reader)  # Read header
            
            for row in reader:
                if row:  # Skip empty rows
                    # Convert to float where possible
                    converted_row = []
                    for val in row:
                        try:
                            converted_row.append(float(val))
                        except ValueError:
                            converted_row.append(val)
                    data['rows'].append(converted_row)
        
        print(f"Loaded {filename}: {len(data['rows'])} data points")
        return data
    
    def get_column(self, data, col_name):
        """Extract column from data"""
        if data is None:
            return None
        
        try:
            col_idx = data['headers'].index(col_name)
            return [row[col_idx] for row in data['rows']]
        except ValueError:
            return None
    
    def plot_resistance_vs_total(self, data, resistance_type, save_file=None):
        """Plot resistance value vs total resistance"""
        if data is None or not data['rows']:
            return
        
        fig, axes = plt.subplots(2, 1, figsize=(12, 10))
        
        # Determine which resistance column to use
        if resistance_type == 'Rx':
            r_col = 'Rx'
            title_suffix = 'Horizontal Resistance (Rx)'
        elif resistance_type == 'Ry':
            r_col = 'Ry'
            title_suffix = 'Vertical Resistance (Ry)'
        elif resistance_type == 'Rz':
            r_col = 'Rz'
            title_suffix = 'Inter-layer Resistance (Rz)'
        else:
            r_col = 'Rx'
            title_suffix = 'All Resistances'
        
        r_values = np.array(self.get_column(data, r_col))
        totalR_values = np.array(self.get_column(data, 'TotalResistance'))
        
        if r_values is None or totalR_values is None:
            print(f"Warning: Cannot plot {resistance_type}, missing data")
            return
        
        # Plot 1: R value vs Total Resistance
        ax1 = axes[0]
        ax1.plot(r_values, totalR_values, 'o-', linewidth=2, markersize=6, 
                color='blue', label='Total Resistance')
        ax1.set_xlabel(f'{title_suffix} (Ω)', fontsize=12, weight='bold')
        ax1.set_ylabel('Total Network Resistance (Ω)', fontsize=12, weight='bold')
        ax1.set_title(f'Total Resistance vs {title_suffix}', fontsize=14, weight='bold')
        ax1.grid(True, alpha=0.3)
        ax1.legend(fontsize=11)
        
        # Add trend line
        z = np.polyfit(r_values, totalR_values, 1)
        p = np.poly1d(z)
        ax1.plot(r_values, p(r_values), "r--", alpha=0.5, linewidth=2, 
                label=f'Linear Fit (slope={z[0]:.2f})')
        ax1.legend(fontsize=11)
        
        # Calculate and display sensitivity
        if len(r_values) >= 2:
            r_change = r_values[-1] - r_values[0]
            totalR_change = totalR_values[-1] - totalR_values[0]
            sensitivity = totalR_change / r_change if r_change != 0 else 0
            
            textstr = f'Sensitivity: {sensitivity:.2f}\n(ΔTotalR / Δ{resistance_type})'
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            ax1.text(0.05, 0.95, textstr, transform=ax1.transAxes, fontsize=11,
                    verticalalignment='top', bbox=props)
        
        # Plot 2: Percentage change
        ax2 = axes[1]
        base_r = r_values[len(r_values)//2]
        base_totalR = totalR_values[len(totalR_values)//2]
        
        r_percent_change = ((r_values - base_r) / base_r) * 100
        totalR_percent_change = ((totalR_values - base_totalR) / base_totalR) * 100
        
        ax2.plot(r_percent_change, totalR_percent_change, 's-', linewidth=2, 
                markersize=6, color='green', label='Total Resistance Change')
        ax2.axhline(y=0, color='k', linestyle='--', alpha=0.3)
        ax2.axvline(x=0, color='k', linestyle='--', alpha=0.3)
        ax2.set_xlabel(f'{resistance_type} Change (%)', fontsize=12, weight='bold')
        ax2.set_ylabel('Total Resistance Change (%)', fontsize=12, weight='bold')
        ax2.set_title(f'Percentage Change: {resistance_type} vs Total Resistance', 
                     fontsize=14, weight='bold')
        ax2.grid(True, alpha=0.3)
        ax2.legend(fontsize=11)
        
        max_change = max(abs(r_percent_change.max()), abs(r_percent_change.min()))
        ax2.plot([-max_change, max_change], [-max_change, max_change], 
                'r--', alpha=0.3, label='1:1 Reference')
        ax2.legend(fontsize=11)
        
        plt.tight_layout()
        
        if save_file:
            os.makedirs(os.path.dirname(save_file) if os.path.dirname(save_file) else '.', exist_ok=True)
            plt.savefig(save_file, dpi=300, bbox_inches='tight')
            print(f"Saved: {save_file}")
        else:
            plt.show()
        plt.close()
    
    def plot_comparison(self, save_file=None):
        """Compare sensitivity of Rx, Ry, Rz"""
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        
        sensitivities = {}
        labels = []
        values = []
        
        for r_type in ['Rx', 'Ry', 'Rz']:
            filename = f'sensitivity_{r_type.lower()}.csv'
            data = self.load_csv(filename)
            if data is not None and len(data['rows']) >= 2:
                if r_type == 'Rx':
                    r_col = 'Rx'
                elif r_type == 'Ry':
                    r_col = 'Ry'
                else:
                    r_col = 'Rz'
                
                r_values = np.array(self.get_column(data, r_col))
                totalR_values = np.array(self.get_column(data, 'TotalResistance'))
                
                if r_values is not None and totalR_values is not None and len(r_values) >= 2:
                    r_change = r_values[-1] - r_values[0]
                    totalR_change = totalR_values[-1] - totalR_values[0]
                    sensitivity = totalR_change / r_change if r_change != 0 else 0
                    
                    sensitivities[r_type] = sensitivity
                    labels.append(r_type)
                    values.append(sensitivity)
        
        if not values:
            print("No data available for comparison")
            return
        
        # Bar chart
        colors = ['blue', 'green', 'orange']
        bars = ax.bar(labels, values, color=colors, alpha=0.7, edgecolor='black', linewidth=2)
        
        # Add value labels on bars
        for bar, val in zip(bars, values):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{val:.2f}', ha='center', va='bottom', fontsize=12, weight='bold')
        
        ax.set_ylabel('Sensitivity (ΔTotalR / ΔR)', fontsize=12, weight='bold')
        ax.set_title('Sensitivity Comparison: Rx vs Ry vs Rz', fontsize=14, weight='bold')
        ax.grid(True, alpha=0.3, axis='y')
        
        max_sens = max(values)
        min_sens = min(values)
        textstr = f'Highest: {labels[values.index(max_sens)]} ({max_sens:.2f})\n'
        textstr += f'Lowest: {labels[values.index(min_sens)]} ({min_sens:.2f})'
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=11,
               verticalalignment='top', bbox=props)
        
        plt.tight_layout()
        
        if save_file:
            os.makedirs(os.path.dirname(save_file) if os.path.dirname(save_file) else '.', exist_ok=True)
            plt.savefig(save_file, dpi=300, bbox_inches='tight')
            print(f"Saved: {save_file}")
        else:
            plt.show()
        plt.close()
    
    def analyze_and_explain(self, data, resistance_type):
        """Analyze results and generate explanation"""
        if data is None or not data['rows']:
            return "No data available for analysis."
        
        explanation = []
        explanation.append(f"\n{'='*60}")
        explanation.append(f"Analysis: {resistance_type} Sensitivity")
        explanation.append(f"{'='*60}\n")
        
        # Determine R column
        if resistance_type == 'Rx':
            r_col = 'Rx'
        elif resistance_type == 'Ry':
            r_col = 'Ry'
        elif resistance_type == 'Rz':
            r_col = 'Rz'
        else:
            r_col = 'Rx'
        
        r_values = np.array(self.get_column(data, r_col))
        totalR_values = np.array(self.get_column(data, 'TotalResistance'))
        
        if r_values is None or totalR_values is None or len(r_values) < 2:
            return "Insufficient data for analysis."
        
        # Calculate sensitivity
        r_initial = r_values[0]
        r_final = r_values[-1]
        r_change = r_final - r_initial
        r_change_percent = (r_change / r_initial) * 100 if r_initial != 0 else 0
        
        totalR_initial = totalR_values[0]
        totalR_final = totalR_values[-1]
        totalR_change = totalR_final - totalR_initial
        totalR_change_percent = (totalR_change / totalR_initial) * 100 if totalR_initial != 0 else 0
        
        sensitivity = totalR_change / r_change if r_change != 0 else 0
        
        explanation.append(f"1. Parameter Range:")
        explanation.append(f"   {resistance_type}: {r_initial:.6f} Ω → {r_final:.6f} Ω")
        explanation.append(f"   Change: {r_change:.6f} Ω ({r_change_percent:+.2f}%)")
        explanation.append("")
        
        explanation.append(f"2. Total Resistance Response:")
        explanation.append(f"   TotalR: {totalR_initial:.6f} Ω → {totalR_final:.6f} Ω")
        explanation.append(f"   Change: {totalR_change:.6f} Ω ({totalR_change_percent:+.2f}%)")
        explanation.append("")
        
        explanation.append(f"3. Sensitivity:")
        explanation.append(f"   Sensitivity = ΔTotalR / Δ{resistance_type}")
        explanation.append(f"              = {totalR_change:.6f} / {r_change:.6f}")
        explanation.append(f"              = {sensitivity:.6f}")
        explanation.append("")
        
        # Interpretation
        explanation.append("4. Interpretation:")
        if sensitivity > 10:
            explanation.append(f"   ⚠️  HIGH SENSITIVITY ({sensitivity:.2f})")
            explanation.append("   → This R value has STRONG impact on total resistance")
            explanation.append("   → Accurate R value modeling is CRITICAL")
            explanation.append("   → Small errors in R will cause large errors in results")
        elif sensitivity > 5:
            explanation.append(f"   ⚡ MODERATE SENSITIVITY ({sensitivity:.2f})")
            explanation.append("   → This R value has MODERATE impact on total resistance")
            explanation.append("   → Accurate R value modeling is IMPORTANT")
            explanation.append("   → Some tolerance for R value errors is acceptable")
        else:
            explanation.append(f"   ✓ LOW SENSITIVITY ({sensitivity:.2f})")
            explanation.append("   → This R value has WEAK impact on total resistance")
            explanation.append("   → Simplified R value modeling may be acceptable")
            explanation.append("   → Larger tolerance for R value errors is acceptable")
        explanation.append("")
        
        # Relationship analysis
        explanation.append("5. Relationship Analysis:")
        correlation = np.corrcoef(r_values, totalR_values)[0, 1]
        explanation.append(f"   Correlation: {correlation:.4f}")
        if abs(correlation) > 0.95:
            explanation.append("   → Strong linear relationship")
        elif abs(correlation) > 0.8:
            explanation.append("   → Moderate linear relationship")
        else:
            explanation.append("   → Weak or non-linear relationship")
        explanation.append("")
        
        return "\n".join(explanation)
    
    def generate_report(self):
        """Generate comprehensive analysis report"""
        report = []
        report.append("="*80)
        report.append("PDN SENSITIVITY ANALYSIS REPORT")
        report.append("="*80)
        report.append("")
        
        # Analyze each resistance type
        for r_type in ['Rx', 'Ry', 'Rz', 'All']:
            filename = f'sensitivity_{r_type.lower()}.csv'
            data = self.load_csv(filename)
            
            if data is not None:
                report.append(self.analyze_and_explain(data, r_type))
                
                # Generate plots
                if r_type != 'All':
                    self.plot_resistance_vs_total(data, r_type, 
                                                  f'{self.output_dir}/plot_{r_type.lower()}_vs_total.png')
        
        # Comparison plot
        self.plot_comparison(f'{self.output_dir}/plot_sensitivity_comparison.png')
        
        # Summary
        report.append("\n" + "="*80)
        report.append("SUMMARY AND RECOMMENDATIONS")
        report.append("="*80)
        report.append("")
        report.append("Based on the sensitivity analysis:")
        report.append("")
        report.append("1. Identify which R value (Rx, Ry, Rz) has the highest sensitivity")
        report.append("2. For high sensitivity R values:")
        report.append("   - Invest more resources in accurate modeling")
        report.append("   - Use detailed extraction methods")
        report.append("   - Validate R values carefully")
        report.append("")
        report.append("3. For low sensitivity R values:")
        report.append("   - Simplified models may be acceptable")
        report.append("   - Some tolerance for errors is acceptable")
        report.append("   - Focus resources on high-sensitivity parameters")
        report.append("")
        
        return "\n".join(report)


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Plot and analyze sensitivity analysis results')
    parser.add_argument('--output-dir', type=str, default='plots', 
                       help='Directory to save plots (default: plots)')
    parser.add_argument('--report', type=str, default='sensitivity_report.txt',
                       help='Output report filename (default: sensitivity_report.txt)')
    parser.add_argument('--no-plots', action='store_true',
                       help='Skip generating plots, only generate report')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Create plotter
    plotter = SensitivityPlotter(output_dir=args.output_dir)
    
    # Check if CSV files exist
    csv_files = ['sensitivity_rx.csv', 'sensitivity_ry.csv', 'sensitivity_rz.csv', 'sensitivity_all.csv']
    found_files = [f for f in csv_files if os.path.exists(f)]
    
    if not found_files:
        print("Error: No CSV files found!")
        print("Please run the sensitivity analysis first:")
        print("  ./pdn_sensitivity")
        sys.exit(1)
    
    print(f"Found {len(found_files)} CSV file(s)")
    print("="*60)
    
    # Generate plots and report
    if not args.no_plots:
        print("Generating plots...")
    
    report = plotter.generate_report()
    
    # Save report
    with open(args.report, 'w') as f:
        f.write(report)
    
    print("\n" + "="*60)
    print("Analysis Complete!")
    print("="*60)
    print(f"Report saved to: {args.report}")
    if not args.no_plots:
        print(f"Plots saved to: {args.output_dir}/")
    print("\nReport Preview:")
    print(report[:1000] + "..." if len(report) > 1000 else report)


if __name__ == '__main__':
    main()

