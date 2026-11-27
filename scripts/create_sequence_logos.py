#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Analyze sequence features from CSV files and generate sequence logos.

This script scans the input directory for files ending with _length.csv,
extracts the first column (sequences), builds a sequence logo using logomaker,
and saves results as PDF files.
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
    Parse command-line arguments
    '''
    parser = argparse.ArgumentParser(description='Analyze sequence features and generate sequence logos')
    parser.add_argument('-i', '--input_dir', required=True, help='Input directory containing files ending with _length.csv')
    parser.add_argument('-l', '--length_filter', type=int, default=None, help='Optional: analyze only sequences of a specific length')
    parser.add_argument('-n', '--top_n', type=int, default=100, help='Number of sequences to use for logos (default: first 100)')
    parser.add_argument('-p', '--position', type=str, default='all', choices=['all', '5p', '3p'], help="Region to analyze: 'all', "
                        "'5p' (5' end), or '3p' (3' end)")
    parser.add_argument('--logfc_filter', type=str, default=None, help='Optional: filter sequences by logFC (e.g., ">2" or "<0")')
    return parser.parse_args()

def create_position_frequency_matrix(sequences, position='all'):
    '''
    Create a position frequency matrix
    
    Args:
        sequences: list of sequences
        position: region to analyze: 'all' (full), '5p' (5' end), or '3p' (3' end)
    
    Returns:
        A position frequency matrix (DataFrame)
    '''
    # Determine sequence region
    if position == 'all':
        # Use full sequences
        pass
    elif position == '5p':
        # Use first 10 bases if sequence is long enough
        sequences = [seq[:10] if len(seq) >= 10 else seq for seq in sequences]
    elif position == '3p':
        # Use last 10 bases if sequence is long enough
        sequences = [seq[-10:] if len(seq) >= 10 else seq for seq in sequences]
    
    # Find maximum sequence length
    max_length = max(len(seq) for seq in sequences)
    
    # Initialize counts matrix
    counts_matrix = {}
    for base in 'ACGTU':
        counts_matrix[base] = [0] * max_length
    
    # Count base occurrences at each position
    for seq in sequences:
        for i, base in enumerate(seq):
            base = base.upper()
            if base in counts_matrix:
                counts_matrix[base][i] += 1
    
    # Convert to DataFrame
    counts_df = pd.DataFrame(counts_matrix)
    
    # Merge 'U' into 'T' if present
    if 'U' in counts_df.columns:
        if 'T' in counts_df.columns:
            counts_df['T'] = counts_df['T'] + counts_df['U']
        else:
            counts_df['T'] = counts_df['U']
        counts_df = counts_df.drop('U', axis=1)
    
    # Totals per position
    totals = counts_df.sum(axis=1)
    
    # Convert counts to frequencies
    freq_matrix = counts_df.div(totals, axis=0).fillna(0)
    
    return freq_matrix

def create_sequence_logo(freq_matrix, output_file, title=None):
    '''
    Create sequence logo with logomaker
    
    Args:
        freq_matrix: position frequency matrix
        output_file: path to output PDF
        title: plot title
    '''
    # Font settings (include common CJK-safe fonts for broader support)
    plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans', 'Bitstream Vera Sans', 'sans-serif']
    plt.rcParams['axes.unicode_minus'] = False  # Ensure minus sign renders properly
    
    plt.figure(figsize=(10, 3))
    
    # Create logo object (use default system font)
    logo = logomaker.Logo(freq_matrix, shade_below=.5, fade_below=.5)
    
    # Styling
    logo.style_spines(visible=False)
    logo.style_xticks(rotation=0, fmt='%d', anchor=0)
    
    # Y-axis limits and labels
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
    
    print(f"Generated sequence logo: {output_file}")

def process_csv_file(file_path, output_dir, length_filter=None, top_n=100, position='all', logfc_filter=None):
    '''
    Process a CSV file and generate a sequence logo
    
    Args:
        file_path: path to CSV
        output_dir: output directory
        length_filter: filter by sequence length
        top_n: number of sequences to use for logo
        position: region to analyze
        logfc_filter: logFC filter string (e.g., ">2")
    '''
    print(f"\nProcessing file: {os.path.basename(file_path)}")
    try:
        # Read CSV file
        df = pd.read_csv(file_path)
        
        # Validate CSV format
        if len(df.columns) < 1:
            print(f"Warning: {file_path} requires at least 1 column")
            return
        
        # First column contains sequences
        sequences = df.iloc[:, 0].tolist()
        
        # If a logFC column exists, filter by logFC
        if logfc_filter is not None and 'logFC' in df.columns:
            match = re.match(r'([<>])([\d.-]+)', str(logfc_filter))
            if match:
                operator, value = match.groups()
                value = float(value)
                if operator == '>':
                    df = df[df['logFC'] > value]
                elif operator == '<':
                    df = df[df['logFC'] < value]
                sequences = df.iloc[:, 0].tolist()
        
        # If a Length column exists, filter by length
        if length_filter is not None and 'Length' in df.columns:
            df = df[df['Length'] == length_filter]
            sequences = df.iloc[:, 0].tolist()
        
        # Ensure we have sequences after filtering
        if len(sequences) == 0:
            print(f"Warning: {file_path} has no sequences after applying filters")
            return
        
        # Use only the first N sequences
        sequences = sequences[:min(top_n, len(sequences))]
        
        # Build position frequency matrix
        freq_matrix = create_position_frequency_matrix(sequences, position)
        
        # Generate output filename
        base_name = os.path.basename(file_path).replace('_length.csv', '')
        position_suffix = '' if position == 'all' else f'_{position}'
        length_suffix = '' if length_filter is None else f'_length{length_filter}'
        logfc_suffix = '' if logfc_filter is None else f'_logFC{logfc_filter}'
        output_file = os.path.join(output_dir, f"{base_name}{position_suffix}{length_suffix}{logfc_suffix}_logo.pdf")
        
        # Build plot title
        title_parts = []
        title_parts.append(base_name)
        if position != 'all':
            end_label = "5' end" if position == '5p' else ("3' end" if position == '3p' else position)
            title_parts.append(end_label)
        if length_filter is not None:
            title_parts.append(f"Length={length_filter}nt")
        if logfc_filter is not None:
            title_parts.append(f"logFC{logfc_filter}")
        title = ' '.join(title_parts)
        
        # Create sequence logo
        create_sequence_logo(freq_matrix, output_file, title)
        
    except Exception as e:
        print(f"Error while processing {file_path}: {str(e)}")

def main():
    args = parse_args()
    
    # Get all files ending with _length.csv in the input directory
    input_pattern = os.path.join(args.input_dir, '*_length.csv')
    csv_files = glob.glob(input_pattern)
    
    if not csv_files:
        print(f"Error: no files ending with _length.csv found in {args.input_dir}")
        sys.exit(1)
    
    # Ensure output directory exists
    os.makedirs(args.input_dir, exist_ok=True)
    
    # Process each CSV file
    for file_path in csv_files:
        process_csv_file(
            file_path, 
            args.input_dir, 
            args.length_filter, 
            args.top_n, 
            args.position,
            args.logfc_filter
        )
    
    print(f"Completed processing {len(csv_files)} file(s)")

if __name__ == '__main__':
    main()
