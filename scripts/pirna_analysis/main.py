#!/usr/bin/env python3
"""
Main program module - handles CLI argument parsing and overall workflow control.
"""

import os
import sys
import time
import argparse
import multiprocessing as mp
from tqdm import tqdm
import pandas as pd

# Import project modules
from pirna_analysis.data_parser import parse_map_file_optimized
from pirna_analysis.data_analyzer import (
    analyze_length_distribution,
    analyze_5prime_shift,
    analyze_3prime_position,
    analyze_unique_reads_length,
    analyze_base_composition_by_position
)
from pirna_analysis.visualization import (
    plot_combined_distributions,
    plot_combined_unique_reads_length,
    plot_base_composition_logo
)
from pirna_analysis.utils import extract_sample_name

def process_file(map_file, ref_file, output_prefix, verbose=False):
    """
    Process a single map file (supports multiprocessing).
    Only performs data processing and saving; no plotting here.
    """
    sample_name = extract_sample_name(map_file)
    if verbose:
        print(f"Processing sample {sample_name}: {map_file}")
    
    # Use the optimized parsing function
    read_info, read_to_pirna_count, read_to_pirnas, filter_stats = parse_map_file_optimized(map_file, ref_file)
    
    # Base composition analysis
    df_base = analyze_base_composition_by_position(read_info, read_to_pirna_count, read_to_pirnas)
    
    if verbose:
        total_reads, filtered_reads = filter_stats
        print(f"  Found {len(read_info)} unique reads in {sample_name}")
        print(f"  Filtered {filtered_reads}/{total_reads} reads with start position > 9 (adjusted threshold > 6)")
    
    # Analyze length distribution
    df_length = analyze_length_distribution(read_info, read_to_pirna_count)
    
    # Analyze 5' end shift distribution
    df_5prime = analyze_5prime_shift(read_info, read_to_pirna_count, read_to_pirnas)
    
    # Analyze 3' end relative position distribution
    df_3prime = analyze_3prime_position(read_info, read_to_pirna_count, read_to_pirnas)
    
    # Analyze unique reads length distribution
    df_unique = analyze_unique_reads_length(read_info)
    
    # Ensure output directory exists
    os.makedirs(output_prefix, exist_ok=True)
    
    # Save results
    df_length.to_csv(f"{output_prefix}/{sample_name}_length_distribution.csv")
    df_5prime.to_csv(f"{output_prefix}/{sample_name}_5prime_shift.csv")
    df_3prime.to_csv(f"{output_prefix}/{sample_name}_3prime_position.csv")
    df_unique.to_csv(f"{output_prefix}/{sample_name}_unique_reads_length.csv")
    df_base.to_csv(f"{output_prefix}/{sample_name}_base_composition.csv")
    
    # Return sample name
    return sample_name

def main():
    # Configure CLI arguments
    parser = argparse.ArgumentParser(description="Analyze piRNA length distribution, 5' end shift, 3' end relative position, and unique reads length distribution")
    parser.add_argument('-i', '--input', nargs='+', required=True, help='One or more input .map files')
    parser.add_argument('-o', '--output', required=True, help='Output directory or prefix')
    parser.add_argument('-v', '--verbose', action='store_true', help='Show verbose processing information')
    parser.add_argument('-p', '--prefix', default='', help='Output filename prefix (default: empty)')
    parser.add_argument('-r', '--ref', default='/home/gyk/project/ld_pirna/data/bmo.v3.0.fa', help='Reference piRNA FASTA file path')
    parser.add_argument('--threads', type=int, default=0, help='Number of worker processes (0 uses all available CPU cores)')
    
    # Parse CLI arguments
    args = parser.parse_args()
    
    # Collect input files and output prefix
    map_files = args.input
    output_prefix = args.output
    verbose = args.verbose
    file_prefix = args.prefix
    
    # Determine number of processes
    threads = args.threads
    if threads <= 0:
        threads = mp.cpu_count()
    
    # Create output directory if needed
    output_dir = os.path.dirname(output_prefix) if os.path.dirname(output_prefix) else output_prefix
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    print(f"Using {threads} CPU core(s) for processing")
    
    # Process multiple samples
    all_data = {}
    all_unique_data = {}
    start_time = time.time()
    
    # Parallel processing path
    if len(map_files) > 1 and threads > 1:
        print(f"Processing {len(map_files)} files in parallel...")
        # Create process pool
        pool = mp.Pool(min(threads, len(map_files)))
        # Prepare arguments for worker function
        process_args = [(map_file, args.ref, output_prefix, verbose) for map_file in map_files]
        sample_names = pool.starmap(process_file, process_args)
        pool.close()
        pool.join()
    else:
        sample_names = []
        for map_file in map_files:
            sample_name = process_file(map_file, args.ref, output_prefix, verbose)
            sample_names.append(sample_name)

    # Aggregate results: read generated CSVs
    for sample_name in sample_names:
        sample_prefix = os.path.join(args.output, sample_name)
        df_length = pd.read_csv(f"{sample_prefix}_length_distribution.csv")
        df_5prime = pd.read_csv(f"{sample_prefix}_5prime_shift.csv")
        df_3prime = pd.read_csv(f"{sample_prefix}_3prime_position.csv")
        df_unique = pd.read_csv(f"{sample_prefix}_unique_reads_length.csv")
        df_base = pd.read_csv(f"{sample_prefix}_base_composition.csv", index_col=0)
        all_data[sample_name] = (df_length, df_5prime, df_3prime)
        all_unique_data[sample_name] = df_unique
        # Plot base composition logo
        plot_base_composition_logo(df_base, f"{sample_prefix}_base_logo.pdf", title=f"{sample_name} base composition")

    # Always plot combined distributions, even for a single sample
    if verbose:
        print("Generating combined distribution plots...")
    
    # Plot combined distributions; export each subplot as a separate figure
    plot_combined_distributions(all_data, f"{output_prefix}/{file_prefix}")
    
    # Plot combined unique reads length distribution
    plot_combined_unique_reads_length(all_unique_data, f"{output_prefix}/{file_prefix}")
    
    print(f"Analysis completed. Results saved to: {output_prefix}")
    print(f"Total time: {time.time() - start_time:.2f} seconds")

if __name__ == "__main__":
    main()
