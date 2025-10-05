#!/usr/bin/env python3
"""
Visualization module - generate various plots.
"""

import matplotlib.pyplot as plt
# Removed unused imports (ticker, numpy)
import os
import logomaker

def plot_combined_distributions(all_data, output_prefix):
    """Plot distributions for multiple samples and export three figures separately"""
    # Use a light theme
    plt.style.use('seaborn-v0_8-pastel')

    # Light color palette and styles - extended to support more samples
    colors = ['#8dd3c7', '#fb8072', '#bebada', '#80b1d3', '#fdb462', '#b3de69', '#fccde5', '#d9d9d9', '#bc80bd', '#ccebc5']
    linestyles = ['-', '--', '-.', ':', '-', '--', '-.', ':', '-', '--']
    markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h']

    # Compute optimal bar width
    num_samples = len(all_data)
    bar_width = 0.8 / num_samples if num_samples > 0 else 0.25

    # Simplify sample names by trimming suffixes
    simplified_sample_names = {}
    for sample_name in all_data.keys():
        # Remove suffixes like ".fa.collapsed"
        simple_name = sample_name
        for suffix in ['.map', '.fa.collapsed.map', '.fa']:
            if simple_name.endswith(suffix):
                simple_name = simple_name[:-len(suffix)]
                break
        simplified_sample_names[sample_name] = simple_name

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_prefix) if os.path.dirname(output_prefix) else output_prefix, exist_ok=True)

    # 1. Length distribution (y-axis as percentage)
    plt.figure(figsize=(10, 6))
    for i, (sample_name, (df_length, _, _)) in enumerate(all_data.items()):
        total_count = df_length['Count'].sum()
        percent = (df_length['Count'] / total_count) * 100 if total_count > 0 else df_length['Count']
        x = df_length['Length'] + i * bar_width - (num_samples * bar_width / 2) + (bar_width / 2)
        plt.bar(x, percent, width=bar_width, color=colors[i % len(colors)],
                label=simplified_sample_names[sample_name])

    # Ensure integer ticks on x-axis
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(1))

    plt.xlabel('Length (nt)')
    plt.ylabel('Percentage of reads (%)')
    plt.title('Length Distribution')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_prefix}length_distribution.pdf", format='pdf')
    plt.close()

    # 2. 5' end shift distribution
    plt.figure(figsize=(10, 6))
    for i, (sample_name, (_, df_shift, _)) in enumerate(all_data.items()):
        # Keep only -2, -1, 0, 1, 2
        filtered_df = df_shift[df_shift['Shift'].isin([-2, -1, 0, 1, 2])].copy()
        plt.plot(filtered_df['Shift'], filtered_df['Percentage'],
                 marker=markers[i % len(markers)], linestyle=linestyles[i % len(linestyles)],
                 color=colors[i % len(colors)], linewidth=2,
                 label=simplified_sample_names[sample_name])

    plt.xlabel("5' position shift")
    plt.ylabel('% of analyzed piRNA loci')
    plt.title("5' End Shift")
    # Show ticks only at -2, -1, 0, 1, 2
    plt.xticks([-2, -1, 0, 1, 2])
    # Ensure y-axis starts at 0
    plt.ylim(bottom=0)
    # Ensure integer ticks on x-axis
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(1))
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_prefix}5prime_shift.pdf", format='pdf')
    plt.close()

    # 3. 3' end relative position distribution
    plt.figure(figsize=(10, 6))
    for i, (sample_name, (_, _, df_position)) in enumerate(all_data.items()):
        plt.plot(df_position['Position'], df_position['Percentage'],
                 marker=markers[i % len(markers)], linestyle=linestyles[i % len(linestyles)],
                 color=colors[i % len(colors)], linewidth=2,
                 label=simplified_sample_names[sample_name])

    plt.xlabel("3' relative position")
    plt.ylabel('% of analyzed piRNA loci')
    plt.title("3' End Position")
    # Ensure that the x-axis scale is an integer.
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(1))
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_prefix}3prime_position.pdf", format='pdf')
    plt.close()

def plot_base_composition_logo(df_base, output_file, title=None):
    """
    Plot base composition sequence logo using logomaker.
    df_base: rows are positions; columns are A/T/C/G/U (percentages 0-100 or fractions)
    output_file: output image path (pdf or png)
    title: optional title
    """
    # Ensure 'Position' is index; keep only A/C/G/U/T columns
    if 'Position' in df_base.columns:
        df_plot = df_base.set_index('Position')
    else:
        df_plot = df_base.copy()
    base_cols = [col for col in df_plot.columns if col in ('A','C','G','U','T')]
    df_plot = df_plot[base_cols]
    # Convert percentage to fraction if needed
    if df_plot.values.max() > 1.1:
        df_plot = df_plot / 100
    plt.figure(figsize=(min(16, 0.4*len(df_plot)), 4))
    logomaker.Logo(df_plot, shade_below=.5, fade_below=.5, font_name='Arial')
    plt.xlabel('Position')
    plt.ylabel('Frequency')
    if title:
        plt.title(title)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def plot_combined_unique_reads_length(all_unique_data, output_prefix):
    """Plot combined unique-read length distributions across samples"""
    # Use a light theme
    plt.style.use('seaborn-v0_8-pastel')

    # Create figure
    plt.figure(figsize=(12, 7))

    # Light color palette and styles - extended to support more samples
    colors = ['#8dd3c7', '#fb8072', '#bebada', '#80b1d3', '#fdb462', '#b3de69', '#fccde5', '#d9d9d9', '#bc80bd', '#ccebc5']

    # Compute optimal bar width
    num_samples = len(all_unique_data)
    bar_width = 0.8 / num_samples if num_samples > 0 else 0.25

    # Plot length distribution (percentage)
    for i, (sample_name, df_unique_length) in enumerate(all_unique_data.items()):
        # Simplify sample names by trimming suffixes
        simplified_name = sample_name
        # Remove common suffixes
        for suffix in ['.map', '.fa.rmdup.map', '.fa.collapsed.map', '.fa.collapsed.no-dust.map', '.collapsed.no-dust.map', '.no-dust.map', '.fa.collapsed.no-dust', '.collapsed.no-dust', '.no-dust', '.fa']:
            if simplified_name.endswith(suffix):
                simplified_name = simplified_name[:-len(suffix)]
                break

        # Remove sequencing-related suffixes
        patterns_to_remove = ['_1.fq', '_2.fq', '.fq', '_1.fastq', '_2.fastq', '.fastq']
        for pattern in patterns_to_remove:
            if pattern in simplified_name:
                simplified_name = simplified_name.replace(pattern, '')

        # If name contains "-KD", keep only that segment
        if '-KD' in simplified_name:
            parts = simplified_name.split('-KD')
            if len(parts) > 1:
                simplified_name = parts[0] + '-KD' + parts[1].split('_')[0]

        total_unique = df_unique_length['UniqueCount'].sum()
        percent_unique = (df_unique_length['UniqueCount'] / total_unique) * 100 if total_unique > 0 else df_unique_length['UniqueCount']
        x = df_unique_length['Length'] + i * bar_width - (num_samples * bar_width / 2) + (bar_width / 2)
        plt.bar(x, percent_unique, width=bar_width, color=colors[i % len(colors)], label=simplified_name)

    plt.xlabel('Length (nt)')
    plt.ylabel('Percentage of unique sequences (%)')
    plt.title('Unique Reads Length Distribution Comparison')
    # 确保x轴刻度为整数
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(1))
    plt.legend()

    # Tight layout
    plt.tight_layout()

    # Save figure as PDF with concise filename
    plt.savefig(f"{output_prefix}unique_length_distribution.pdf", format='pdf')
    plt.close()
