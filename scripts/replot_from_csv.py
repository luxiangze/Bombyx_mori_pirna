#!/usr/bin/env python3
"""
Replot piRNA analysis figures from CSV files to avoid re-parsing map files.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import glob


# Use a light color theme consistent with the original script
plt.style.use('seaborn-v0_8-pastel')

def extract_sample_name(file_path):
    """Extract sample name from file path"""
    file_name = os.path.basename(file_path)
    # Remove common suffixes
    for suffix in ['_length_distribution.csv', '_5prime_shift.csv', '_3prime_position.csv', '_unique_reads_length.csv']:
        if file_name.endswith(suffix):
            file_name = file_name[:-len(suffix)]
            break
    # Remove other common suffixes
    for suffix in ['.fa.collapsed', '.map', '.no-dust']:
        if file_name.endswith(suffix):
            file_name = file_name[:-len(suffix)]
    return file_name

def load_csv_data(input_dir):
    """Load all CSV files under the given directory"""
    data = {}

    # Find CSV files for each plot type
    length_files = glob.glob(os.path.join(input_dir, '*_length_distribution.csv'))
    shift_files = glob.glob(os.path.join(input_dir, '*_5prime_shift.csv'))
    position_files = glob.glob(os.path.join(input_dir, '*_3prime_position.csv'))
    unique_length_files = glob.glob(os.path.join(input_dir, '*_unique_reads_length.csv'))

    # Load data for each sample
    for length_file in length_files:
        sample_name = extract_sample_name(length_file)

        # Find the corresponding CSV files for the sample
        shift_file = None
        for f in shift_files:
            if extract_sample_name(f) == sample_name:
                shift_file = f
                break

        position_file = None
        for f in position_files:
            if extract_sample_name(f) == sample_name:
                position_file = f
                break

        unique_length_file = None
        for f in unique_length_files:
            if extract_sample_name(f) == sample_name:
                unique_length_file = f
                break

        # Load dataframes
        df_length = pd.read_csv(length_file)
        df_shift = pd.read_csv(shift_file) if shift_file else None
        df_position = pd.read_csv(position_file) if position_file else None
        df_unique_length = pd.read_csv(unique_length_file) if unique_length_file else None

        data[sample_name] = (df_length, df_shift, df_position, df_unique_length)

    return data

def plot_combined_distributions(all_data, output_prefix):
    """Plot combined distributions"""
    # Use a light color theme consistent with the original script
    plt.style.use('seaborn-v0_8-pastel')

    # Light color palettes and styles - extended to support more samples
    colors = ['#8dd3c7', '#fb8072', '#bebada', '#80b1d3', '#fdb462', '#b3de69', '#fccde5', '#d9d9d9', '#bc80bd', '#ccebc5']
    markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h']
    linestyles = ['-', '--', '-.', ':', '-', '--', '-.', ':', '-', '--']

    # Compute the number of samples and bar width
    num_samples = len(all_data)
    bar_width = 0.8 / num_samples if num_samples > 1 else 0.6

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

    # 1. Length distribution
    plt.figure(figsize=(10, 6))
    for i, (sample_name, (df_length, _, _, _)) in enumerate(all_data.items()):
        total_count = df_length['Count'].sum()
        percentage = df_length['Count'] / total_count * 100 if total_count > 0 else df_length['Count']
        x = df_length['Length'] + i * bar_width - (num_samples * bar_width / 2) + (bar_width / 2)
        plt.bar(x, percentage, width=bar_width, color=colors[i % len(colors)],
                label=simplified_sample_names[sample_name])

    # Ensure integer ticks on x-axis
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(1))

    plt.xlabel('Length (nt)')
    plt.ylabel('Percentage (%)')
    plt.title('Length Distribution')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_prefix}length_distribution.pdf", format='pdf')
    plt.close()

    # 2. 5' end shift distribution
    plt.figure(figsize=(10, 6))
    for i, (sample_name, (_, df_shift, _, _)) in enumerate(all_data.items()):
        if df_shift is None:
            continue
        # Keep only five values: -2, -1, 0, 1, 2
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
    # 确保x轴刻度为整数
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(1))
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_prefix}5prime_shift.pdf", format='pdf')
    plt.close()

    # 3. 3' end relative position distribution
    plt.figure(figsize=(10, 6))
    for i, (sample_name, (_, _, df_position, _)) in enumerate(all_data.items()):
        if df_position is None:
            continue
        plt.plot(df_position['Position'], df_position['Percentage'],
                 marker=markers[i % len(markers)], linestyle=linestyles[i % len(linestyles)],
                 color=colors[i % len(colors)], linewidth=2,
                 label=simplified_sample_names[sample_name])

    plt.xlabel("3' relative position")
    plt.ylabel('% of analyzed piRNA loci')
    plt.title("3' End Position")
    # 确保x轴刻度为整数
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(1))
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_prefix}3prime_position.pdf", format='pdf')
    plt.close()

    # 4. Unique reads length distribution
    plt.figure(figsize=(10, 6))
    for i, (sample_name, (_, _, _, df_unique_length)) in enumerate(all_data.items()):
        if df_unique_length is None:
            continue

        # Simplify sample name by trimming suffixes
        simplified_name = simplified_sample_names[sample_name]
        # Remove sequencing-related suffixes
        patterns_to_remove = ['_1.fq', '_2.fq', '.fq', '_1.fastq', '_2.fastq', '.fastq']
        for pattern in patterns_to_remove:
            if pattern in simplified_name:
                simplified_name = simplified_name.replace(pattern, '')

        # If sample name contains "-KD", keep only that part
        if '-KD' in simplified_name:
            parts = simplified_name.split('-KD')
            if len(parts) > 1:
                simplified_name = parts[0] + '-KD' + parts[1].split('_')[0]

        total_unique = df_unique_length['UniqueCount'].sum()
        percentage = df_unique_length['UniqueCount'] / total_unique * 100 if total_unique > 0 else df_unique_length['UniqueCount']
        x = df_unique_length['Length'] + i * bar_width - (num_samples * bar_width / 2) + (bar_width / 2)
        plt.bar(x, percentage, width=bar_width, color=colors[i % len(colors)],
                label=simplified_name)

    plt.xlabel('Length (nt)')
    plt.ylabel('Percentage (%)')
    plt.title('Unique Reads Length Distribution')
    # 确保x轴刻度为整数
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(1))
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_prefix}unique_length_distribution.pdf", format='pdf')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Replot piRNA analysis figures from CSV files')
    parser.add_argument('-i', '--input', required=True, help='Input directory containing CSV files')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('-p', '--prefix', default='', help='Output filename prefix')
    args = parser.parse_args()

    # Ensure output directory exists
    os.makedirs(args.output, exist_ok=True)

    # Load CSV data
    all_data = load_csv_data(args.input)

    if not all_data:
        print("No valid CSV files found!")
        return

    # Set output prefix
    output_prefix = os.path.join(args.output, args.prefix)

    # Plot combined distributions
    plot_combined_distributions(all_data, output_prefix)

    print(f"Successfully replotted figures; saved to directory: {args.output}")

if __name__ == "__main__":
    main()
