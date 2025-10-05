#!/usr/bin/env python3
"""
Data parsing module - responsible for reading and parsing map files and reference piRNA files.
"""

import os
import time
import pandas as pd
from collections import defaultdict
from tqdm import tqdm
from Bio import SeqIO

def parse_map_file_optimized(map_file, ref_pirna_file, chunk_size=10000):
    """
    Optimized version of parse_map_file using chunked processing and efficient data structures.
    
    Args:
        map_file: path to map file
        ref_pirna_file: path to reference piRNA FASTA
        chunk_size: number of rows per chunk to process
    
    Returns:
        read_info, read_to_pirna_count, read_to_pirnas, filter_stats
    """
    # Store number of reference piRNAs each read maps to
    read_to_pirna_count = defaultdict(int)
    # Store count and length per read sequence
    read_info = {}
    # Store mapped reference piRNAs and positions per read sequence
    read_to_pirnas = defaultdict(list)
    
    # Track filtered sequences
    filtered_count = 0
    total_count = 0
    
    # Timing
    start_time = time.time()
    
    print(f"Start parsing reference piRNA file: {ref_pirna_file}")
    # Read sequences from reference piRNA file
    pirna_sequences = {}
    try:
        # Use iterator mode to reduce memory usage
        for record in SeqIO.parse(ref_pirna_file, "fasta"):
            pirna_sequences[record.id] = str(record.seq)
        print(f"Loaded {len(pirna_sequences)} reference piRNA sequences in {time.time() - start_time:.2f} s")
    except Exception as e:
        print(f"Warning: failed to read reference piRNA file: {e}")
        print("Will use target_seq as a fallback")
    
    # Find corresponding .fa.collapsed file to get sequence counts
    collapsed_file = map_file.replace('.map', '')
    if not os.path.exists(collapsed_file):
        # Try alternative filename patterns
        possible_collapsed_files = [
            map_file.replace('.fa.rmdup.map', '.fa.collapsed'),
            map_file.replace('.map', '.fa.collapsed'),
            os.path.join(os.path.dirname(map_file), os.path.basename(map_file).split('.')[0] + '.fa.collapsed')
        ]
        
        for possible_file in possible_collapsed_files:
            if os.path.exists(possible_file):
                collapsed_file = possible_file
                print(f"Found collapsed file: {collapsed_file}")
                break
    
    # Read sequence counts and sequences from collapsed FASTA via Biopython
    seq_counts = {}
    read_sequences = {}
    if os.path.exists(collapsed_file):
        collapsed_start_time = time.time()
        print(f"Parsing collapsed file: {collapsed_file}")
        try:
            # Count lines to initialize progress bar
            with open(collapsed_file, 'r') as f:
                line_count = sum(1 for _ in f)
            
            record_count = 0
            for record in tqdm(SeqIO.parse(collapsed_file, "fasta"), total=line_count//2, desc="Reading sequences"):
                # Extract count from FASTA header
                count = 1
                if '-' in record.id:
                    try:
                        count = int(record.id.split('-')[-1])
                    except ValueError:
                        pass
                
                # Store sequence and its count
                seq = str(record.seq)
                seq_counts[seq] = count
                
                # Map read ID to sequence string
                read_sequences[record.id] = seq
                record_count += 1
            
            print(f"Loaded {len(seq_counts)} sequences with counts in {time.time() - collapsed_start_time:.2f} s")
            print(f"Stored {len(read_sequences)} read ID to sequence mappings")
        except Exception as e:
            print(f"Warning: failed to parse collapsed file: {e}")
    else:
        print(f"Warning: collapsed file not found: {collapsed_file}")
    
    # Read map file
    map_start_time = time.time()
    print(f"Start parsing map file: {map_file}")
    try:
        # Read TSV map file in chunks with pandas to reduce memory usage
        # First, get column count
        header = pd.read_csv(map_file, sep='\t', nrows=1, header=None)
        column_count = len(header.columns)
        
        # Set column names
        column_names = ['piRNA_id', 'trans_coord', 'reads_id']
        
        # Count lines to configure progress bar
        with open(map_file, 'r') as f:
            line_count = sum(1 for _ in f) - 1  # minus header line
        
        # Read and process by chunks
        reader = pd.read_csv(map_file, sep='\t', skiprows=1, header=None, names=column_names, chunksize=chunk_size)
        
        # Progress bar
        with tqdm(total=line_count, desc="Processing map file") as pbar:
            for chunk in reader:
                chunk_size = len(chunk)
                total_count += chunk_size
                pbar.update(chunk_size)
                
                # Process each row
                for _, row in chunk.iterrows():
                    piRNA_id = row['piRNA_id']  # reference piRNA ID
                    trans_coord = int(row['trans_coord'])  # offset position
                    reads_id = row['reads_id']  # read sequence ID
                    
                    # Fetch read sequence and its count
                    read_seq = read_sequences.get(reads_id, "")
                    if not read_seq:
                        # Skip if read sequence not found
                        continue
                    count = seq_counts.get(read_seq, 1)
                    
                    # Store read info
                    read_key = read_seq
                    read_info[read_key] = (count, len(read_seq))
                    
                    # Update number of mapped reference piRNAs for this read
                    read_to_pirna_count[read_key] += 1
                    
                    # Append mapped reference piRNA and position
                    read_to_pirnas[read_key].append((piRNA_id, trans_coord))
        
        print(f"Parsed {total_count} alignment records in {time.time() - map_start_time:.2f} s")
        print(f"Filtered {filtered_count} sequences with start position > 6")
    except Exception as e:
        print(f"Warning: failed to parse map file: {e}")
    
    # Return parsed results and filter statistics
    filter_stats = (total_count, filtered_count)
    return read_info, read_to_pirna_count, read_to_pirnas, filter_stats
