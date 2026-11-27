#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Convert a CSV file to FASTA format.

Assumes the first three columns are: sequence, control reads, treatment reads.
Generates two FASTA files. Each record ID is in the form "index-reads".
Output filenames are derived from the two sample names in the CSV filename.
'''

import os
import sys
import csv
import re
import argparse

def parse_args():
    '''
    Parse command-line arguments
    '''
    parser = argparse.ArgumentParser(description='Convert CSV to FASTA format')
    parser.add_argument('csv_file', help='Path to input CSV file')
    parser.add_argument('-o', '--output_dir', 
                      default='/home/gyk/project/ld_pirna/data/bm_smRNA_pirna_diffrencial', 
                      help='Output directory (default: /home/gyk/project/ld_pirna/data/bm_smRNA_pirna_diffrencial)')
    return parser.parse_args()

def extract_sample_names(filename):
    '''
    Extract sample names from CSV filename.
    Example: from "Control_vs_SUGP1-KD2_edger_results_with_counts_FDR_0.05_length.csv"
    extract "Control" and "SUGP1-KD2".
    '''
    base_name = os.path.basename(filename)
    match = re.search(r'(.+?)_vs_(.+?)_', base_name)
    if match:
        return match.group(1), match.group(2)
    else:
        # Fallback if names cannot be extracted
        return 'sample1', 'sample2'

def csv_to_fasta(csv_file, output_dir='/home/gyk/project/ld_pirna/data/bm_smRNA_pirna_diffrencial'):
    '''
    Convert CSV to FASTA files (control and treatment)
    '''
    # Extract sample names from filename
    sample1, sample2 = extract_sample_names(csv_file)
    
    # Extract the name prefix, e.g., "Control_vs_SUGP1-KD2"
    base_name = os.path.basename(csv_file)
    match = re.search(r'(.+?)_edger', base_name)
    if match:
        name_prefix = match.group(1)
    else:
        # Fallback to sample names
        name_prefix = f"{sample1}_vs_{sample2}"
    
    # Create output directory if it does not exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Define output file paths with concise names
    fasta_file1 = os.path.join(output_dir, f"{name_prefix}-control.fa.collapsed")
    fasta_file2 = os.path.join(output_dir, f"{name_prefix}-treatment.fa.collapsed")
    
    # Open output files
    with open(fasta_file1, 'w') as f1, open(fasta_file2, 'w') as f2:
        # Read CSV file
        with open(csv_file, 'r') as csvfile:
            reader = csv.reader(csvfile)
            # Skip header line
            header = next(reader)
            
            # Validate CSV format
            if len(header) < 3:
                print(f"Error: CSV format invalid; expected at least 3 columns but got {len(header)}", file=sys.stderr)
                sys.exit(1)
            
            # Iterate through rows
            for i, row in enumerate(reader, 1):
                if len(row) < 3:
                    print(f"Warning: row {i} is incomplete; skipped", file=sys.stderr)
                    continue
                
                # Extract sequence and read counts
                sequence = row[0]
                control_reads = int(float(row[1])) if row[1] else 0
                treatment_reads = int(float(row[2])) if row[2] else 0
                
                # Write FASTA records
                if control_reads > 0:
                    f1.write(f">{i}-{control_reads}\n{sequence}\n")
                if treatment_reads > 0:
                    f2.write(f">{i}-{treatment_reads}\n{sequence}\n")
    
    print(f"Generated FASTA files:\n{fasta_file1}\n{fasta_file2}")
    return fasta_file1, fasta_file2

def main():
    args = parse_args()
    csv_to_fasta(args.csv_file, args.output_dir)

if __name__ == '__main__':
    main()
