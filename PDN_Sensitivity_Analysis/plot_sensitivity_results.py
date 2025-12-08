#!/usr/bin/env python3
"""
Plot and Analyze Sensitivity Analysis Results

This script reads CSV files from sensitivity analysis and creates visualizations
with detailed explanations of the phenomena observed.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

# Set style for better-looking plots
plt.style.use('seaborn-v0_8-darkgrid' if 'seaborn-v0_8-darkgrid' in plt.style.available else 'default')
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 10

class SensitivityPlotter:
    def __init__(self, output_dir='.'):
        self.output_dir = output_dir
        self.results = {}
        
    def load_csv(self, filename):
        """Load CSV file and return DataFrame"""
        if not os.path.exists(filename):
            print(f"Warning: {filename} not found, skipping...")
            return None
        
        df = pd.read_csv(filename)
        print(f"Loaded {filename}: {len(df)} data points")
        return df
    
    def plot_resistance_vs_total(self, df, resistance_type, save_file=None):
        """
        Plot resistance value vs total resistance
        Shows how individual R value affects total network resistance
        """
        if df is None or df.empty:
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
            r_col = 'Rx'  # Default
            title_suffix = 'All Resistances'
        
        # Plot 1: R value vs Total Resistance
        ax1 = axes[0]
        ax1.plot(df[r_col], df['TotalResistance'], 'o-', linewidth=2, markersize=6, 
                color='blue', label='Total Resistance')
        ax1.set_xlabel(f'{title_suffix} (Ω)', fontsize=12, weight='bold')
        ax1.set_ylabel('Total Network Resistance (Ω)', fontsize=12, weight='bold')
        ax1.set_title(f'Total Resistance vs {title_suffix}', fontsize=14, weight='bold')
        ax1.grid(True, alpha=0.3)
        ax1.legend(fontsize=11)
        
        # Add trend line
        z = np.polyfit(df[r_col], df['TotalResistance'], 1)
        p = np.poly1d(z)
        ax1.plot(df[r_col], p(df[r_col]), "r--", alpha=0.5, linewidth=2, 
                label=f'Linear Fit (slope={z[0]:.2f})')
        ax1.legend(fontsize=11)
        
        # Calculate and display sensitivity
        if len(df) >= 2:
            r_change = df[r_col].iloc[-1] - df[r_col].iloc[0]
            totalR_change = df['TotalResistance'].iloc[-1] - df['TotalResistance'].iloc[0]
            sensitivity = totalR_change / r_change if r_change != 0 else 0
            
            # Add text box with sensitivity info
            textstr = f'Sensitivity: {sensitivity:.2f}\n(ΔTotalR / Δ{resistance_type})'
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            ax1.text(0.05, 0.95, textstr, transform=ax1.transAxes, fontsize=11,
                    verticalalignment='top', bbox=props)
        
        # Plot 2: Percentage change
        ax2 = axes[1]
        base_r = df[r_col].iloc[len(df)//2]  # Use middle value as base
        base_totalR = df['TotalResistance'].iloc[len(df)//2]
        
        r_percent_change = ((df[r_col] - base_r) / base_r) * 100
        totalR_percent_change = ((df['TotalResistance'] - base_totalR) / base_totalR) * 100
        
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
        
        # Add 1:1 reference line
        max_change = max(abs(r_percent_change.max()), abs(r_percent_change.min()))
        ax2.plot([-max_change, max_change], [-max_change, max_change], 
                'r--', alpha=0.3, label='1:1 Reference')
        ax2.legend(fontsize=11)
        
        plt.tight_layout()
        
        if save_file:
            plt.savefig(save_file, dpi=300, bbox_inches='tight')
            print(f"Saved: {save_file}")
        else:
            plt.show()
        plt.close()
    
    def plot_branch_statistics(self, df, resistance_type, save_file=None):
        """Plot branch resistance statistics"""
        if df is None or df.empty:
            return
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        if resistance_type == 'Rx':
            r_col = 'Rx'
        elif resistance_type == 'Ry':
            r_col = 'Ry'
        elif resistance_type == 'Rz':
            r_col = 'Rz'
        else:
            r_col = 'Rx'
        
        # Plot 1: Average branch resistance
        ax1 = axes[0, 0]
        ax1.plot(df[r_col], df['AvgBranchResistance'], 'o-', linewidth=2, 
                markersize=5, color='blue')
        ax1.set_xlabel(f'{resistance_type} (Ω)', fontsize=11)
        ax1.set_ylabel('Average Branch Resistance (Ω)', fontsize=11)
        ax1.set_title('Average Branch Resistance', fontsize=12, weight='bold')
        ax1.grid(True, alpha=0.3)
        
        # Plot 2: Min branch resistance
        ax2 = axes[0, 1]
        ax2.plot(df[r_col], df['MinBranchResistance'], 's-', linewidth=2, 
                markersize=5, color='green')
        ax2.set_xlabel(f'{resistance_type} (Ω)', fontsize=11)
        ax2.set_ylabel('Min Branch Resistance (Ω)', fontsize=11)
        ax2.set_title('Minimum Branch Resistance', fontsize=12, weight='bold')
        ax2.grid(True, alpha=0.3)
        
        # Plot 3: Max branch resistance
        ax3 = axes[1, 0]
        ax3.plot(df[r_col], df['MaxBranchResistance'], '^-', linewidth=2, 
                markersize=5, color='red')
        ax3.set_xlabel(f'{resistance_type} (Ω)', fontsize=11)
        ax3.set_ylabel('Max Branch Resistance (Ω)', fontsize=11)
        ax3.set_title('Maximum Branch Resistance', fontsize=12, weight='bold')
        ax3.grid(True, alpha=0.3)
        
        # Plot 4: Resistance range (max - min)
        ax4 = axes[1, 1]
        resistance_range = df['MaxBranchResistance'] - df['MinBranchResistance']
        ax4.plot(df[r_col], resistance_range, 'd-', linewidth=2, 
                markersize=5, color='orange')
        ax4.set_xlabel(f'{resistance_type} (Ω)', fontsize=11)
        ax4.set_ylabel('Resistance Range (Ω)', fontsize=11)
        ax4.set_title('Branch Resistance Range (Max - Min)', fontsize=12, weight='bold')
        ax4.grid(True, alpha=0.3)
        
        plt.suptitle(f'Branch Resistance Statistics: {resistance_type} Analysis', 
                    fontsize=14, weight='bold', y=1.02)
        plt.tight_layout()
        
        if save_file:
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
            df = self.load_csv(filename)
            if df is not None and len(df) >= 2:
                if r_type == 'Rx':
                    r_col = 'Rx'
                elif r_type == 'Ry':
                    r_col = 'Ry'
                else:
                    r_col = 'Rz'
                
                r_change = df[r_col].iloc[-1] - df[r_col].iloc[0]
                totalR_change = df['TotalResistance'].iloc[-1] - df['TotalResistance'].iloc[0]
                sensitivity = totalR_change / r_change if r_change != 0 else 0
                
                sensitivities[r_type] = sensitivity
                labels.append(r_type)
                values.append(sensitivity)
        
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
        
        # Add interpretation
        max_sens = max(values) if values else 0
        min_sens = min(values) if values else 0
        
        if max_sens > 0:
            textstr = f'Highest: {labels[values.index(max_sens)]} ({max_sens:.2f})\n'
            textstr += f'Lowest: {labels[values.index(min_sens)]} ({min_sens:.2f})'
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=11,
                   verticalalignment='top', bbox=props)
        
        plt.tight_layout()
        
        if save_file:
            plt.savefig(save_file, dpi=300, bbox_inches='tight')
            print(f"Saved: {save_file}")
        else:
            plt.show()
        plt.close()
    
    def analyze_and_explain(self, df, resistance_type):
        """Analyze results and generate explanation"""
        if df is None or df.empty:
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
        
        # Calculate sensitivity
        if len(df) >= 2:
            r_initial = df[r_col].iloc[0]
            r_final = df[r_col].iloc[-1]
            r_change = r_final - r_initial
            r_change_percent = (r_change / r_initial) * 100 if r_initial != 0 else 0
            
            totalR_initial = df['TotalResistance'].iloc[0]
            totalR_final = df['TotalResistance'].iloc[-1]
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
            correlation = df[r_col].corr(df['TotalResistance'])
            explanation.append(f"   Correlation: {correlation:.4f}")
            if abs(correlation) > 0.95:
                explanation.append("   → Strong linear relationship")
            elif abs(correlation) > 0.8:
                explanation.append("   → Moderate linear relationship")
            else:
                explanation.append("   → Weak or non-linear relationship")
            explanation.append("")
            
            # Impact assessment
            explanation.append("6. Impact Assessment:")
            if r_change_percent != 0:
                impact_ratio = totalR_change_percent / abs(r_change_percent)
                explanation.append(f"   Impact Ratio: {impact_ratio:.4f}")
                explanation.append(f"   (TotalR changes {impact_ratio:.2f}x the rate of {resistance_type})")
                if impact_ratio > 1.0:
                    explanation.append("   → TotalR is MORE sensitive than the R value itself")
                elif impact_ratio < 1.0:
                    explanation.append("   → TotalR is LESS sensitive than the R value itself")
                else:
                    explanation.append("   → TotalR changes proportionally with R value")
            explanation.append("")
        
        # Branch statistics
        if 'AvgBranchResistance' in df.columns:
            explanation.append("7. Branch Resistance Statistics:")
            explanation.append(f"   Average: {df['AvgBranchResistance'].mean():.6f} Ω")
            explanation.append(f"   Min: {df['MinBranchResistance'].min():.6f} Ω")
            explanation.append(f"   Max: {df['MaxBranchResistance'].max():.6f} Ω")
            explanation.append(f"   Range: {df['MaxBranchResistance'].max() - df['MinBranchResistance'].min():.6f} Ω")
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
            df = self.load_csv(filename)
            
            if df is not None:
                report.append(self.analyze_and_explain(df, r_type))
                
                # Generate plots
                if r_type != 'All':
                    self.plot_resistance_vs_total(df, r_type, 
                                                  f'{self.output_dir}/plot_{r_type.lower()}_vs_total.png')
                    self.plot_branch_statistics(df, r_type, 
                                               f'{self.output_dir}/plot_{r_type.lower()}_branch_stats.png')
        
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
        report.append("4. Overall:")
        report.append("   - High sensitivity → Finding 'good R values' is CRITICAL")
        report.append("   - Low sensitivity → Finding 'good R values' is less critical")
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


