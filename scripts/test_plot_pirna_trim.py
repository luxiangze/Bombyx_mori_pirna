#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
import matplotlib.pyplot as plt

# Import original plot_pirna_trim and wrap a fixed version
from pirna_analysis.visualization import plot_pirna_trim as original_plot_pirna_trim  # noqa: F401

# Fixed version of plot_pirna_trim
def fixed_plot_pirna_trim(all_pirna_trim_data, output_prefix):
    """Fixed plot_pirna_trim that corrects colorbar issues"""
    # Use a light theme
    plt.style.use('seaborn-v0_8-pastel')
    
    # Gather samples
    sample_names = list(all_pirna_trim_data.keys())
    n_samples = len(sample_names)
    
    # Require at least 2 samples for comparison
    if n_samples < 2:
        print("Warning: at least 2 samples are required for comparison")
        return
    
    # Simplify sample names
    simplified_sample_names = {}
    for sample_name in sample_names:
        # Remove common suffixes
        simple_name = sample_name
        for suffix in ['.map', '.fa.rmdup.map', '.fa.collapsed.map', '.fa.collapsed.no-dust.map', 
                      '.collapsed.no-dust.map', '.no-dust.map', '.fa.collapsed.no-dust', 
                      '.collapsed.no-dust', '.no-dust', '.fa']:
            if simple_name.endswith(suffix):
                simple_name = simple_name[:-len(suffix)]
                break
        
        # Remove sequencing-related suffixes
        patterns_to_remove = ['_1.fq', '_2.fq', '.fq', '_1.fastq', '_2.fastq', '.fastq']
        for pattern in patterns_to_remove:
            if pattern in simple_name:
                simple_name = simple_name.replace(pattern, '')
        
        # If contains "-KD", keep that segment
        if '-KD' in simple_name:
            parts = simple_name.split('-KD')
            if len(parts) > 1:
                simple_name = parts[0] + '-KD' + parts[1].split('_')[0]
        
        simplified_sample_names[sample_name] = simple_name
    
    # Colormap for fold change
    cmap = plt.cm.coolwarm
    
    # Plot pairwise comparisons and save each as a separate file
    for i in range(n_samples):
        for j in range(i+1, n_samples):
            # Get simplified names
            sample1_name = simplified_sample_names[sample_names[i]]
            sample2_name = simplified_sample_names[sample_names[j]]
            
            # Create figure/axes
            fig, ax = plt.subplots(figsize=(6, 6), dpi=100)
            
            # Fetch data
            df1 = all_pirna_trim_data[sample_names[i]]
            df2 = all_pirna_trim_data[sample_names[j]]
            
            # Merge by piRNA_ID
            merged_df = pd.merge(df1, df2, on='piRNA_ID', how='inner', suffixes=('_1', '_2'))
            
            # Compute color by fold change (add epsilon to avoid division by zero)
            epsilon = 1e-10
            fold_changes = merged_df['Trim_index_2'] / (merged_df['Trim_index_1'] + epsilon)
            
            # Determine fc range for symmetric mapping around 1
            fc_min = fold_changes.min()
            fc_max = fold_changes.max()
            
            # Ensure symmetry around 1
            if fc_min < 1 and 1/fc_min > fc_max:
                vmax = 1/fc_min
                vmin = fc_min
            else:
                vmax = fc_max
                vmin = 1/fc_max
            
            # Clamp range to avoid extremes
            vmax = min(vmax, 10)
            vmin = max(vmin, 0.1)
            
            # Normalizer
            norm = plt.Normalize(vmin=vmin, vmax=vmax)
            
            # Map each point to a color
            colors = [cmap(norm(fc)) for fc in fold_changes]
            
            # Scatter plot
            ax.scatter(merged_df['Trim_index_1'], merged_df['Trim_index_2'], 
                      c=colors, alpha=0.7, s=15, edgecolors='none')
            
            # Diagonal y=x
            max_val = max(merged_df['Trim_index_1'].max(), merged_df['Trim_index_2'].max())
            min_val = min(merged_df['Trim_index_1'].min(), merged_df['Trim_index_2'].min())
            # Ensure min=1 for log scale
            min_val = max(1, min_val)
            
            # Draw diagonal
            ax.plot([min_val, max_val], [min_val, max_val], 'k-', linewidth=1)
            
            # Log scales
            ax.set_xscale('log')
            ax.set_yscale('log')
            
            # Limits
            ax.set_xlim(min_val, max_val*1.1)
            ax.set_ylim(min_val, max_val*1.1)
            
            # Axis labels
            ax.set_xlabel(f'{sample1_name}\n(Normalized reads)')
            ax.set_ylabel(f'{sample2_name}\n(Normalized reads)')
            
            # ScalarMappable for colorbar
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            
            # Colorbar on the figure
            cbar = fig.colorbar(sm, ax=ax)
            cbar.set_label(f'Fold change ({vmin:.2f}-{vmax:.2f})')
            
            # Layout
            plt.tight_layout()
            
            # Save figure
            output_filename = f"{output_prefix}{sample1_name}_vs_{sample2_name}"
            plt.savefig(f"{output_filename}.pdf", format='pdf', bbox_inches='tight')
            plt.close()

def main():
    """Test plot_pirna_trim"""
    # Set input/output dirs
    input_dir = "/home/gyk/project/ld_pirna/results/pirna_analysis_20250405_164942"
    output_dir = "/home/gyk/project/ld_pirna/results/test_plot_pirna_trim"
    
    # Ensure output dir exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Read all pirna_trim data
    all_pirna_trim_data = {}
    sample_files = [f for f in os.listdir(input_dir) if f.endswith("_pirna_trim.csv")]
    
    print("Found the following sample files:")
    for sample_file in sample_files:
        sample_name = sample_file.replace("_pirna_trim.csv", "")
        print(f"  - {sample_file}")
        
        # Read CSV
        file_path = os.path.join(input_dir, sample_file)
        df = pd.read_csv(file_path)
        
        # Basic info
        print(f"    Rows: {len(df)}")
        print(f"    Columns: {', '.join(df.columns)}")
        
        # Store
        all_pirna_trim_data[sample_name] = df
    
    print(f"\nTotal samples found: {len(all_pirna_trim_data)}")
    
    # Call fixed version
    print("\nStart plotting piRNA trimming comparison...")
    fixed_plot_pirna_trim(all_pirna_trim_data, f"{output_dir}/")
    
    print(f"\nPlotting completed. Results saved to: {output_dir}")

if __name__ == "__main__":
    main()
