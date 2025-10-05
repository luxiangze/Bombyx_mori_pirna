#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse
import os

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Plot sequence length distribution from a TSV file.')
    parser.add_argument('input_file', help='Path to the length distribution TSV file')
    parser.add_argument('-o', '--output', help='Output directory for the plots', default='.')
    parser.add_argument('--dpi', type=int, default=300, help='DPI for output images')
    parser.add_argument('--figsize', type=str, default='12,8', help='Figure size in inches, format: width,height')
    parser.add_argument('--log', action='store_true', help='Use log scale for y-axis')
    parser.add_argument('--min-length', type=int, default=18, help='Minimum sequence length to include')
    parser.add_argument('--max-length', type=int, default=50, help='Maximum sequence length to include')
    return parser.parse_args()

def read_data(input_file, min_length=18, max_length=50):
    """Read the length distribution data from the TSV file."""
    # Read the TSV file
    df = pd.read_csv(input_file, sep='\t')
    
    # Filter by length range
    df = df[(df['Length'] >= min_length) & (df['Length'] <= max_length)]
    
    # Set Length as index for easier plotting
    df.set_index('Length', inplace=True)
    
    return df

def convert_to_percentage(df):
    """Convert raw counts to percentages for each sample."""
    # Calculate the sum for each column (total reads per sample)
    column_sums = df.sum()
    
    # Divide each value by the column sum and multiply by 100 to get percentage
    df_percent = df.div(column_sums) * 100
    
    return df_percent

def plot_all_samples(df, output_dir, figsize=(12, 8), dpi=300, log_scale=False):
    """Plot all samples in a single figure."""
    plt.figure(figsize=figsize)
    
    # Set the style
    sns.set_style("whitegrid")
    
    # Plot each sample
    for column in df.columns:
        plt.plot(df.index, df[column], marker='o', linewidth=2, markersize=4, label=column)
    
    # Set plot properties
    plt.title('Sequence Length Distribution', fontsize=16)
    plt.xlabel('Sequence Length (nt)', fontsize=14)
    plt.ylabel('Percentage (%)', fontsize=14)
    
    if log_scale:
        plt.yscale('log')
        plt.ylabel('Percentage (%) - log scale', fontsize=14)
    
    plt.xticks(df.index, rotation=0)
    plt.grid(True, alpha=0.3)
    plt.legend(loc='best', fontsize=10)
    
    # Tight layout
    plt.tight_layout()
    
    # Save the figure
    output_file = os.path.join(output_dir, 'all_samples_length_distribution_percent.png')
    plt.savefig(output_file, dpi=dpi)
    print(f"Saved plot to {output_file}")
    
    plt.close()

def plot_by_group(df, output_dir, figsize=(12, 8), dpi=300, log_scale=False):
    """Group samples by their prefix and plot each group."""
    # Extract sample prefixes (assuming format is PREFIX-REP)
    prefixes = [col.split('-')[0] for col in df.columns]
    unique_prefixes = list(set(prefixes))
    
    # Create a plot for each group
    for prefix in unique_prefixes:
        plt.figure(figsize=figsize)
        
        # Get columns for this prefix
        group_columns = [col for col in df.columns if col.startswith(prefix)]
        
        # Plot each sample in the group
        for column in group_columns:
            plt.plot(df.index, df[column], marker='o', linewidth=2, markersize=4, label=column)
        
        # Set plot properties
        plt.title(f'Sequence Length Distribution - {prefix}', fontsize=16)
        plt.xlabel('Sequence Length (nt)', fontsize=14)
        plt.ylabel('Percentage (%)', fontsize=14)
        
        if log_scale:
            plt.yscale('log')
            plt.ylabel('Percentage (%) - log scale', fontsize=14)
        
        plt.xticks(df.index, rotation=0)
        plt.grid(True, alpha=0.3)
        plt.legend(loc='best', fontsize=10)
        
        # Tight layout
        plt.tight_layout()
        
        # Save the figure
        output_file = os.path.join(output_dir, f'{prefix}_length_distribution_percent.png')
        plt.savefig(output_file, dpi=dpi)
        print(f"Saved plot to {output_file}")
        
        plt.close()

def plot_heatmap(df, output_dir, figsize=(12, 10), dpi=300):
    """Create a heatmap of the length distribution."""
    plt.figure(figsize=figsize)
    
    # Normalize data for better visualization
    df_norm = df.div(df.max())
    
    # Create heatmap
    sns.heatmap(df_norm.T, cmap='viridis', linewidths=0.5, linecolor='white')
    
    plt.title('Normalized Sequence Length Distribution Heatmap', fontsize=16)
    plt.xlabel('Sequence Length (nt)', fontsize=14)
    plt.ylabel('Sample', fontsize=14)
    
    # Tight layout
    plt.tight_layout()
    
    # Save the figure
    output_file = os.path.join(output_dir, 'length_distribution_heatmap.png')
    plt.savefig(output_file, dpi=dpi)
    print(f"Saved heatmap to {output_file}")
    
    plt.close()

def main():
    """Main function."""
    # Parse command line arguments
    args = parse_args()
    
    # Parse figsize
    figsize = tuple(map(float, args.figsize.split(',')))
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)
    
    # Read data
    df = read_data(args.input_file, args.min_length, args.max_length)
    
    # Convert to percentage
    df_percent = convert_to_percentage(df)
    
    # Create plots
    plot_all_samples(df_percent, args.output, figsize, args.dpi, args.log)
    plot_by_group(df_percent, args.output, figsize, args.dpi, args.log)
    plot_heatmap(df, args.output, figsize, args.dpi)  # Keep heatmap with original values
    
    print("All plots generated successfully!")

if __name__ == "__main__":
    main()
