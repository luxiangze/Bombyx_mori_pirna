#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Generate dual sequence logo plots based on sequence data and expression levels in CSV files

This script analyzes sequence data in input files,
generates two sequence logo plots based on sample expression levels in the second and third columns,
and saves the results in PDF format.
'''

import os
import sys
import glob
import pandas as pd
import matplotlib.pyplot as plt
import logomaker
import argparse
import re

def parse_args():
    '''
    Parse command line arguments
    '''
    parser = argparse.ArgumentParser(description='Generate dual sequence logo plots based on expression levels')
    parser.add_argument('-i', '--input_dir', required=True, help='Input directory containing CSV files ending with length.csv')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory path')
    parser.add_argument('-l', '--length_filter', type=int, default=None, help='Optional, analyze only sequences of specific length')
    parser.add_argument('-n', '--top_n', type=int, default=None, help='Number of sequences used to generate the logo, default is all sequences')
    parser.add_argument('-p', '--position', type=str, default='all', choices=['all', '5p', '3p'], help='Analysis position, options: all (entire sequence), 5p (5\' end), or 3p (3\' end)')
    parser.add_argument('--logfc_filter', type=str, default=None, help='Optional, filter sequences by logFC value, e.g., ">2" means sequences with logFC>2')
    parser.add_argument('--control_name', type=str, default=None, help='Control group name (second column), if not provided, will be extracted from filename')
    parser.add_argument('--treatment_name', type=str, default=None, help='Treatment group name (third column), if not provided, will be extracted from filename')
    parser.add_argument('--min_count', type=int, default=1, help='Minimum expression threshold, sequences below this value will be filtered, default is 0')
    parser.add_argument('--sequence_filter', type=str, default=None, help='Optional, exclude sequences containing specific motif, e.g., "AATT"')
    parser.add_argument('--auto_detect', action='store_true', help='Automatically detect control and treatment names from filename')
    return parser.parse_args()

def create_position_frequency_matrix(sequences, weights=None, position='all'):
    '''
    Create position frequency matrix with weighted sequences

    Parameters:
    sequences - list of sequences
    weights - list of weights for each sequence (e.g., expression levels)
    position - analysis position, 'all' for entire sequence, '5p' for 5' end, '3p' for 3' end

    Returns:
    position frequency matrix
    '''
    # If no weights provided, use equal weights
    if weights is None:
        weights = [1] * len(sequences)

    # Ensure sequences and weights have the same length
    if len(sequences) != len(weights):
        raise ValueError("Sequences and weights must have the same length")

    # Determine which part of sequences to use by position
    if position == 'all':
        # Use entire sequence
        pass
    elif position == '5p':
        # Use first 10 nt if available
        sequences = [seq[:10] if len(seq) >= 10 else seq for seq in sequences]
    elif position == '3p':
        # Use last 10 nt if available
        sequences = [seq[-10:] if len(seq) >= 10 else seq for seq in sequences]

    # Find max length among sequences
    max_length = max(len(seq) for seq in sequences) if sequences else 0

    if max_length == 0:
        return pd.DataFrame()

    # Initialize counts matrix
    counts_matrix = {}
    for base in 'ACGTU':
        counts_matrix[base] = [0] * max_length

    # Accumulate weighted base counts per position
    for seq, weight in zip(sequences, weights):
        for i, base in enumerate(seq):
            base = base.upper()
            if base in counts_matrix:
                counts_matrix[base][i] += weight

    # Convert to DataFrame
    counts_df = pd.DataFrame(counts_matrix)

    # Merge 'U' into 'T' column if present
    if 'U' in counts_df.columns:
        if 'T' in counts_df.columns:
            counts_df['T'] = counts_df['T'] + counts_df['U']
        else:
            counts_df['T'] = counts_df['U']
        counts_df = counts_df.drop('U', axis=1)

    # Total counts per position
    totals = counts_df.sum(axis=1)

    # Normalize to frequency per position
    freq_matrix = counts_df.div(totals, axis=0).fillna(0)

    return freq_matrix, counts_df

def create_sequence_logo(freq_matrix, output_file, title=None):
    '''
    Create sequence logo plot using logomaker

    Parameters:
    freq_matrix - position frequency matrix
    output_file - output file path
    title - chart title
    '''
    if freq_matrix.empty:
        print(f"Warning: Not enough sequence data to generate logo plot: {output_file}")
        return

    # Set font
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Bitstream Vera Sans', 'sans-serif']
    plt.rcParams['axes.unicode_minus'] = False  # Fix minus sign display issue

    plt.figure(figsize=(10, 3))

    # Build Logo using frequency matrix directly; no IC transform.
    # Turn off shading/fade for clarity.
    logo = logomaker.Logo(
        freq_matrix,
        shade_below=0,  # No shadows used
        fade_below=0,
    )   # No fade effect

    # Styling
    logo.style_spines(visible=False)
    logo.style_xticks(rotation=0, fmt='%d', anchor=0)

    # Axis labels and limits
    plt.ylim(0, 1)
    plt.ylabel('Frequency')
    plt.xlabel('Position')

    # Title
    if title:
        plt.title(title)

    # Save figure
    plt.tight_layout()
    plt.savefig(output_file, format='pdf', dpi=300)
    plt.close()

    print(f"Sequence logo plot generated: {output_file}")

def extract_sample_names_from_filename(file_path):
    '''
    Extract control and treatment names from filename

    Parameters:
    file_path - path to CSV file

    Returns:
    control_name, treatment_name - extracted names
    '''
    filename = os.path.basename(file_path)
    # Remove extension and common suffix
    base_name = filename.replace('.csv', '').replace('_length', '')

    # Try to extract sample names from filename pattern like "control_vs_treatment"
    if '_vs_' in base_name:
        parts = base_name.split('_vs_')
        if len(parts) >= 2:
            control_name = parts[0]
            treatment_name = parts[1]
            return control_name, treatment_name

    # Fallback to other common patterns
    patterns = [
        r'(\w+)[_-]to[_-](\w+)',  # pattern: control_to_treatment
        r'(\w+)[_-]and[_-](\w+)',  # pattern: control_and_treatment
        r'(\w+)[_-](\w+)'  # simple two-part pattern
    ]

    for pattern in patterns:
        match = re.search(pattern, base_name)
        if match:
            return match.group(1), match.group(2)

    # Default if extraction fails
    return 'Control', 'Treatment'

def process_csv_file(file_path, args):
    '''
    Process CSV file and generate two sequence logo plots

    Parameters:
    file_path - path to CSV file
    args - command line arguments
    '''
    output_dir = args.output_dir
    length_filter = args.length_filter
    top_n = args.top_n
    position = args.position
    logfc_filter = args.logfc_filter
    sequence_filter = args.sequence_filter
    min_count = args.min_count

    # Auto-detect sample names or use provided ones
    if args.auto_detect or (args.control_name is None and args.treatment_name is None):
        control_name, treatment_name = extract_sample_names_from_filename(file_path)
        print(f"Auto-detected sample names: Control={control_name}, Treatment={treatment_name}")
    else:
        control_name = args.control_name if args.control_name is not None else 'Control'
        treatment_name = args.treatment_name if args.treatment_name is not None else 'Treatment'

    print(f"\nProcessing file: {os.path.basename(file_path)}")
    try:
        # Read CSV
        df = pd.read_csv(file_path)

        # Validate format
        if len(df.columns) < 3:
            print(f"Error: {file_path} has incorrect format, at least 3 columns required (sequence, control expression, treatment expression)")
            return

        # Ensure output dir exists
        os.makedirs(output_dir, exist_ok=True)

        # Column names
        seq_col = df.columns[0]
        control_col = df.columns[1]
        treatment_col = df.columns[2]

        # Apply filters
        filtered_df = df.copy()

        # Exclude sequences containing motif
        if sequence_filter is not None:
            filtered_df = filtered_df[~filtered_df[seq_col].str.contains(sequence_filter, case=False)]
            print(f"Excluded {len(df) - len(filtered_df)} sequences containing '{sequence_filter}'")

        # Filter by logFC if present
        if logfc_filter is not None and 'logFC' in filtered_df.columns:
            match = re.match(r'([<>])([\d.-]+)', str(logfc_filter))
            if match:
                operator, value = match.groups()
                value = float(value)
                if operator == '>':
                    filtered_df = filtered_df[filtered_df['logFC'] > value]
                elif operator == '<':
                    filtered_df = filtered_df[filtered_df['logFC'] < value]

        # Filter by length if present
        if length_filter is not None and 'Length' in filtered_df.columns:
            filtered_df = filtered_df[filtered_df['Length'] == length_filter]

        # Ensure sequences remain after filtering
        if len(filtered_df) == 0:
            print(f"Warning: no sequences remain after filters for {file_path}")
            return

        # Base output name
        filename = os.path.basename(file_path)

        # Simplify base name
        if '_edger_results_with_counts_FDR_' in filename:
            # Extract sample name part
            parts = filename.split('_edger_results_with_counts_FDR_')
            if len(parts) >= 2:
                # Keep first part and FDR value
                fdr_parts = parts[1].split('_')
                if len(fdr_parts) >= 2:
                    base_name = f"{parts[0]}_edger_results_with_counts_FDR_{fdr_parts[0]}"
                else:
                    base_name = parts[0]
            else:
                base_name = filename.replace('.csv', '').replace('_length', '')
        else:
            base_name = filename.replace('.csv', '').replace('_length', '')

        # Simplify treatment name
        if '_edger_results_with_counts_FDR_' in treatment_name:
            treatment_name_parts = treatment_name.split('_edger_results_with_counts_FDR_')
            if len(treatment_name_parts) >= 1:
                treatment_name = treatment_name_parts[0]

        position_suffix = '' if position == 'all' else f'_{position}'
        length_suffix = '' if length_filter is None else f'_length{length_filter}'
        logfc_suffix = '' if logfc_filter is None else f'_logFC{logfc_filter}'

        # Process control group data
        print(f"Processing control group ({control_name}) data...")
        # Filter and sort by control expression
        control_df = filtered_df[filtered_df[control_col] >= min_count].sort_values(by=control_col, ascending=False)

        if len(control_df) > 0:
            # Use all or top N
            if top_n is not None:
                control_df = control_df.head(min(top_n, len(control_df)))

            # Get sequences and weights
            control_sequences = control_df[seq_col].tolist()
            control_weights = control_df[control_col].tolist()

            # Create PFM with weighted sequences
            control_freq_matrix, control_counts = create_position_frequency_matrix(control_sequences, control_weights, position)

            # Build output filenames
            control_output_file = os.path.join(output_dir, f"{base_name}_{control_name}{position_suffix}{length_suffix}{logfc_suffix}_logo.pdf")
            control_output_table = os.path.join(output_dir, f"{base_name}_{control_name}{position_suffix}{length_suffix}{logfc_suffix}_counts.csv")

            # Save counts table
            control_counts.to_csv(control_output_table, index=False)

            # Title pieces
            control_title_parts = [f"{base_name} - {control_name}"]
            if position != 'all':
                control_title_parts.append(f"{position} end")
            if length_filter is not None:
                control_title_parts.append(f"length={length_filter}nt")
            if logfc_filter is not None:
                control_title_parts.append(f"logFC{logfc_filter}")
            control_title = ' '.join(control_title_parts)

            # Create logo
            create_sequence_logo(control_freq_matrix, control_output_file, control_title)
        else:
            print(f"Warning: Control group ({control_name}) has no sequences with expression level ≥{min_count}")

        # Process treatment group data
        print(f"Processing treatment group ({treatment_name}) data...")
        # Filter and sort by treatment expression
        treatment_df = filtered_df[filtered_df[treatment_col] >= min_count].sort_values(by=treatment_col, ascending=False)

        if len(treatment_df) > 0:
            # Use all sequences or top N if specified
            if top_n is not None:
                treatment_df = treatment_df.head(min(top_n, len(treatment_df)))

            # Get sequences and their expression levels
            treatment_sequences = treatment_df[seq_col].tolist()
            treatment_weights = treatment_df[treatment_col].tolist()

            # Create position frequency matrix with weighted sequences
            treatment_freq_matrix, treatment_counts = create_position_frequency_matrix(treatment_sequences, treatment_weights, position)

            # Output filenames
            treatment_output_file = os.path.join(output_dir, f"{base_name}_{treatment_name}{position_suffix}{length_suffix}{logfc_suffix}_logo.pdf")
            treatment_output_table = os.path.join(output_dir, f"{base_name}_{treatment_name}{position_suffix}{length_suffix}{logfc_suffix}_counts.csv")

            # save counts
            treatment_counts.to_csv(treatment_output_table, index=False)

            # Title pieces
            treatment_title_parts = [f"{base_name} - {treatment_name}"]
            if position != 'all':
                treatment_title_parts.append(f"{position} end")
            if length_filter is not None:
                treatment_title_parts.append(f"length={length_filter}nt")
            if logfc_filter is not None:
                treatment_title_parts.append(f"logFC{logfc_filter}")
            treatment_title = ' '.join(treatment_title_parts)

            # Create logo
            create_sequence_logo(treatment_freq_matrix, treatment_output_file, treatment_title)
        else:
            print(f"Warning: Treatment group ({treatment_name}) has no sequences with expression level ≥{min_count}")

    except Exception as e:
        print(f"Error processing file {file_path}: {str(e)}")

def main():
    args = parse_args()

    # Validate input dir
    if not os.path.isdir(args.input_dir):
        print(f"Error: Input directory {args.input_dir} does not exist")
        sys.exit(1)

    # Find all *length.csv files under input dir
    input_pattern = os.path.join(args.input_dir, '*length.csv')
    csv_files = glob.glob(input_pattern)

    if not csv_files:
        print(f"Error: No files ending with length.csv found in {args.input_dir}")
        sys.exit(1)

    print(f"Found {len(csv_files)} CSV files to process")

    # 确保输出目录存在
    os.makedirs(args.output_dir, exist_ok=True)

    # Process each CSV
    for file_path in csv_files:
        process_csv_file(file_path, args)

    print(f"Processing complete. Processed {len(csv_files)} files.")

if __name__ == '__main__':
    main()
