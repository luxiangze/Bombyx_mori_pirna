#!/usr/bin/env python3
"""
Analyze piRNA length distribution, 5' end shift distribution, and 3' end relative position.
Usage:
    python analyze_pirna_distribution.py -i <map_file1> [<map_file2> ...] -o <output_directory>
or:
    python analyze_pirna_distribution.py --input <map_file1> [<map_file2> ...] --output <output_directory>
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import time
import multiprocessing as mp
from collections import defaultdict
from Bio import SeqIO
from tqdm import tqdm

def parse_map_file(map_file, ref_pirna_file='/home/gyk/project/ld_pirna/data/bmo.v3.0.fa'):
    """Parse seqmap-generated map file and extract alignments to reference piRNAs.
    Filter out reads with alignment start position > 6 (original trans_coord > 9).

    Map file example (new format):
    trans_id        trans_coord     target_seq      probe_id        probe_seq       num_mismatch
    piR-bmo-1       1               TCCCGTATGGTCTAGTGGCTA   208883-4    TCCCGTATGGTCTAGTGGCTA   0

    Notes:
    1) We align positions 3–23 of each read; therefore start_pos = trans_coord - 3.
    2) ref_seq should be obtained from ref_pirna_file, not from target_seq.
    3) read_seq is in the collapsed file; use probe_id to retrieve it when available.
    """
    # Store count of reads aligned to each reference piRNA
    read_to_pirna_count = defaultdict(int)
    # Store count and length of each read
    read_info = {}  # {read_key: (count, length)}
    # Store reference piRNA and position for each read
    read_to_pirnas = defaultdict(list)
    # Store read sequence ID and sequence mapping
    read_sequences = {}

    # Record filtered sequence count
    filtered_count = 0
    total_count = 0

    print(f"Parsing reference piRNA file: {ref_pirna_file}")
    # Read sequences from reference piRNA file
    pirna_sequences = {}
    try:
        for record in SeqIO.parse(ref_pirna_file, "fasta"):
            pirna_sequences[record.id] = str(record.seq)
        print(f"Loaded {len(pirna_sequences)} reference piRNA sequences")
    except Exception as e:
        print(f"Warning: error reading reference piRNA file: {e}")
        print("Fallback to target_seq")

    # Find corresponding .fa.collapsed file to get sequence count information
    collapsed_file = map_file.replace('.map', '')
    if not os.path.exists(collapsed_file):
        # Try other possible file name patterns
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

    # Use Biopython to parse collapsed file for sequence count information and sequences
    seq_counts = {}
    if os.path.exists(collapsed_file):
        print(f"Parsing collapsed file: {collapsed_file}")
        try:
            for record in SeqIO.parse(collapsed_file, "fasta"):
                # Extract count information from FASTA header
                count = 1
                if '-' in record.id:
                    try:
                        count = int(record.id.split('-')[-1])
                    except ValueError:
                        pass

                # Store sequence and count
                seq = str(record.seq)
                seq_counts[seq] = count

                # Store sequence ID and sequence mapping
                read_sequences[record.id] = seq

            print(f"Loaded {len(seq_counts)} sequences with counts")
            print(f"Stored {len(read_sequences)} read ID-to-sequence mappings")
        except Exception as e:
            print(f"Warning: error parsing collapsed file: {e}")
    else:
        print(f"Warning: collapsed file not found: {collapsed_file}")

    # Parse map file
    print(f"Parsing map file: {map_file}")
    try:
        # Use pandas to read TSV format map file
        df_map = pd.read_csv(map_file, sep='\t', skiprows=1, header=None)
        # Ensure at least 6 columns
        if df_map.shape[1] < 6:
            raise ValueError(f"Map file format error: fewer than 6 columns ({df_map.shape[1]})")

        # Rename columns for easier understanding
        df_map.columns = ['trans_id', 'trans_coord', 'target_seq', 'probe_id', 'probe_seq', 'num_mismatch'] + \
                          [f'col{i+7}' for i in range(df_map.shape[1]-6)]

        # Process each row
        total_count = len(df_map)
        for _, row in df_map.iterrows():
            pirna_id = row['trans_id']  # Use trans_id as pirna_id
            raw_start_pos = int(row['trans_coord'])

            # Adjust start_pos, because using 3-23nt of each sequence for alignment
            start_pos = raw_start_pos - 3

            # Get ref_seq from reference piRNA file
            if pirna_id in pirna_sequences:
                ref_seq = pirna_sequences[pirna_id]
            else:
                # If not found, use target_seq as fallback
                ref_seq = row['target_seq']
                print(f"Warning: reference piRNA ID not found: {pirna_id}; using target_seq as fallback")

            probe_id = row['probe_id']
            read_seq = row['probe_seq']

            # If probe_id in read_sequences, use mapped sequence
            if probe_id in read_sequences:
                read_seq = read_sequences[probe_id]

            _mismatches = int(row['num_mismatch'])

            # Filter out reads with alignment start position > 6 (original position > 9)
            if raw_start_pos > 9:  # Equivalent to adjusted start_pos > 6
                filtered_count += 1
                continue

            # Extract count information from probe_id, format: ID-COUNT
            read_count = 1  # Default count is 1
            if '-' in probe_id:
                try:
                    read_count = int(probe_id.split('-')[-1])
                except ValueError:
                    pass

            # Get sequence count from seq_counts if available
            if read_seq in seq_counts:
                read_count = seq_counts[read_seq]

            # Use read sequence as key, assume all are forward strand
            strand = '+'
            read_key = f"{read_seq}_{strand}"

            # Update read sequence count and length
            read_info[read_key] = (read_count, len(read_seq))

            # Update read sequence alignment to reference piRNA
            read_to_pirna_count[read_key] += 1
            read_to_pirnas[read_key].append((pirna_id, start_pos, len(ref_seq)))
    except Exception as e:
        print(f"Warning: error parsing map file: {e}")

    # Return filtering statistics for detailed mode display
    filter_stats = (total_count, filtered_count)
    print(f"Parsing complete: total {total_count} records, filtered {filtered_count}")
    return read_info, read_to_pirna_count, read_to_pirnas, filter_stats

def analyze_length_distribution(read_info, read_to_pirna_count):
    """Analyze piRNA length distribution"""
    # Store total count for each length
    length_counts = defaultdict(int)

    # Process each read sequence
    for read_key, (count, length) in read_info.items():
        # Get the number of reference piRNAs aligned to this read sequence
        pirna_count = read_to_pirna_count[read_key]
        # Normalized weight: read count / number of aligned reference piRNAs
        normalized_weight = count / pirna_count

        # Update length distribution
        length_counts[length] += normalized_weight

    # Convert to DataFrame
    lengths = sorted(length_counts.keys())
    counts = [length_counts[length_val] for length_val in lengths]

    df_length = pd.DataFrame({
        'Length': lengths,
        'Count': counts
    })

    return df_length

def analyze_5prime_shift(read_info, read_to_pirna_count, read_to_pirnas):
    """Analyze 5' end shift distribution"""
    # Store normalized count for each shift value
    shift_counts = defaultdict(float)

    # Process each read sequence
    for read_key, (count, read_length) in read_info.items():
        # Get the number of reference piRNAs aligned to this read sequence
        pirna_count = read_to_pirna_count[read_key]
        # Normalized weight: read count / number of aligned reference piRNAs
        normalized_weight = count / pirna_count

        # Update normalized count for each shift value
        for pirna_id, start_pos, ref_length in read_to_pirnas[read_key]:
            # 5' end displacement: starting position - 1 (because positions start from 1)
            # In seqmap format, trans_coord already starts from 1
            shift = start_pos - 1
            shift_counts[shift] += normalized_weight

    # Convert to DataFrame
    shifts = sorted(shift_counts.keys())
    counts = [shift_counts[s] for s in shifts]

    # Calculate total count
    total_count = sum(counts)

    # Calculate percentage
    percentages = [count / total_count * 100 for count in counts] if total_count > 0 else [0] * len(counts)

    df_shift = pd.DataFrame({
        'Shift': shifts,
        'Count': counts,
        'Percentage': percentages
    })

    return df_shift

def analyze_3prime_position(read_info, read_to_pirna_count, read_to_pirnas):
    """Analyze 3' end relative position distribution"""
    # Store normalized count for each 3' end relative position
    position_counts = defaultdict(float)

    # Process each read sequence
    for read_key, (count, read_length) in read_info.items():
        # Get the number of reference piRNAs aligned to this read sequence
        pirna_count = read_to_pirna_count[read_key]
        # Normalized weight: read count / number of aligned reference piRNAs
        normalized_weight = count / pirna_count

        # Update normalized count for each 3' end relative position
        for pirna_id, start_pos, ref_length in read_to_pirnas[read_key]:
            # Calculate 3' end position: starting position + read length - 1
            # In seqmap format, trans_coord already starts from 1
            end_pos = start_pos + read_length - 1
            position_counts[end_pos] += normalized_weight

    # Convert to DataFrame
    positions = sorted(position_counts.keys())
    counts = [position_counts[p] for p in positions]

    # Calculate total count
    total_count = sum(counts)

    # Calculate percentage
    percentages = [count / total_count * 100 for count in counts] if total_count > 0 else [0] * len(counts)

    df_position = pd.DataFrame({
        'Position': positions,
        'Count': counts,
        'Percentage': percentages
    })

    return df_position

def plot_distributions(df_length, df_shift, df_position, sample_name, output_prefix):
    """Plot distribution figures"""
    # Set chart style to pastel
    plt.style.use('seaborn-v0_8-pastel')

    # Create a figure with three subplots
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))

    # Plot length distribution (similar to Figure B)
    ax1.bar(df_length['Length'], df_length['Count'] / 1e6, color='#1f77b4', width=0.8)
    ax1.set_xlabel('Length (nt)')
    ax1.set_ylabel('Reads per million (x 10⁶)')
    ax1.set_title(f'{sample_name} - Length Distribution')

    # Plot 5' end shift distribution (similar to Figure C left)
    ax2.plot(df_shift['Shift'], df_shift['Percentage'], 'o-', color='#d62728', linewidth=2)
    ax2.set_xlabel("5' position shift")
    ax2.set_ylabel('% of analyzed piRNA loci')
    ax2.set_title(f"{sample_name} - 5' End Shift")

    # Plot 3' end relative position distribution (similar to Figure C right)
    ax3.plot(df_position['Position'], df_position['Percentage'], 'o-', color='#ff7f0e', linewidth=2)
    ax3.set_xlabel("3' relative position")
    ax3.set_ylabel('% of analyzed piRNA loci')
    ax3.set_title(f"{sample_name} - 3' End Position")

    # Adjust layout
    plt.tight_layout()

    # Save figure as PDF
    plt.savefig(f"{output_prefix}_{sample_name}_distributions.pdf", format='pdf')
    plt.close()

def plot_combined_distributions(all_data, output_prefix):
    """Plot combined distributions for multiple samples; export three figures"""
    # Set chart style to pastel
    plt.style.use('seaborn-v0_8-pastel')

    # Pastel colors and styles - extended to support more samples
    colors = ['#8dd3c7', '#fb8072', '#bebada', '#80b1d3', '#fdb462', '#b3de69', '#fccde5', '#d9d9d9', '#bc80bd', '#ccebc5']#ccebc5']
    linestyles = ['-', '--', '-.', ':', '-', '--', '-.', ':', '-', '--']
    markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h']

    # Calculate optimal bar width
    num_samples = len(all_data)
    bar_width = 0.8 / num_samples if num_samples > 0 else 0.25

    # Simplify sample names, remove suffixes
    simplified_sample_names = {}
    for sample_name in all_data.keys():
        # Remove suffixes such as ".fa.collapsed"
        simple_name = sample_name
        for suffix in ['.map', '.fa.collapsed.map', '.fa']:
            if simple_name.endswith(suffix):
                simple_name = simple_name[:-len(suffix)]
                break
        simplified_sample_names[sample_name] = simple_name

    # Plot length distribution
    plt.figure(figsize=(10, 6))
    for i, (sample_name, (df_length, _, _)) in enumerate(all_data.items()):
        x = df_length['Length'] + i * bar_width - (num_samples * bar_width / 2) + (bar_width / 2)
        plt.bar(x, df_length['Count'] / 1e6, width=bar_width, color=colors[i % len(colors)],
                label=simplified_sample_names[sample_name])

    plt.xlabel('Length (nt)')
    plt.ylabel('Reads per million (x 10⁶)')
    plt.title('Length Distribution')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_prefix}/length_distribution.pdf", format='pdf')
    plt.close()

    # Plot 5' end shift distribution
    plt.figure(figsize=(10, 6))
    for i, (sample_name, (_, df_shift, _)) in enumerate(all_data.items()):
        # Only keep 0, 1, 2 values
        filtered_df = df_shift[df_shift['Shift'].isin([0, 1, 2])].copy()
        plt.plot(filtered_df['Shift'], filtered_df['Percentage'],
                 marker=markers[i % len(markers)], linestyle=linestyles[i % len(linestyles)],
                 color=colors[i % len(colors)], linewidth=2,
                 label=simplified_sample_names[sample_name])

    plt.xlabel("5' position shift")
    plt.ylabel('% of analyzed piRNA loci')
    plt.title("5' End Shift")
    # Set x axis to only show 0, 1, 2
    plt.xticks([0, 1, 2])
    # Ensure y axis starts from 0
    plt.ylim(bottom=0)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_prefix}/5prime_shift.pdf", format='pdf')
    plt.close()

    # Plot 3' end relative position distribution
    plt.figure(figsize=(10, 6))
    for i, (sample_name, (_, _, df_position)) in enumerate(all_data.items()):
        plt.plot(df_position['Position'], df_position['Percentage'],
                 marker=markers[i % len(markers)], linestyle=linestyles[i % len(linestyles)],
                 color=colors[i % len(colors)], linewidth=2,
                 label=simplified_sample_names[sample_name])

    plt.xlabel("3' relative position")
    plt.ylabel('% of analyzed piRNA loci')
    plt.title("3' End Position")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_prefix}/3prime_position.pdf", format='pdf')
    plt.close()

def analyze_unique_reads_length(read_info):
    """Analyze unique reads length distribution"""
    # Use pandas Series to efficiently handle length counts
    lengths = [length for _, (_, length) in read_info.items()]
    length_series = pd.Series(lengths)

    # Calculate the number of unique sequences for each length
    unique_length_counts = length_series.value_counts().sort_index()

    # Convert to DataFrame
    df_unique_length = pd.DataFrame({
        'Length': unique_length_counts.index,
        'UniqueCount': unique_length_counts.values
    })

    return df_unique_length

def plot_unique_reads_length(df_unique_length, sample_name, output_prefix):
    """Plot unique reads length distribution"""
    # Set chart style to pastel
    plt.style.use('seaborn-v0_8-pastel')

    # Create figure
    plt.figure(figsize=(10, 6))

    # Plot length distribution
    plt.bar(df_unique_length['Length'], df_unique_length['UniqueCount'], color='#b3de69', width=0.8)
    plt.xlabel('Length (nt)')
    plt.ylabel('Number of unique sequences')
    plt.title(f'{sample_name} - Unique Reads Length Distribution')

    # Adjust layout
    plt.tight_layout()

    # Save figure as PDF
    plt.savefig(f"{output_prefix}_{sample_name}_unique_length_distribution.pdf", format='pdf')
    plt.close()

def plot_combined_unique_reads_length(all_unique_data, output_prefix):
    """Plot combined unique reads length distribution for multiple samples"""
    # Set chart style to pastel
    plt.style.use('seaborn-v0_8-pastel')

    # Create figure
    plt.figure(figsize=(12, 7))

    # Pastel colors and styles - extended to support more samples
    colors = ['#8dd3c7', '#fb8072', '#bebada', '#80b1d3', '#fdb462', '#b3de69', '#fccde5', '#d9d9d9', '#bc80bd', '#ccebc5']

    # Calculate optimal bar width
    num_samples = len(all_unique_data)
    bar_width = 0.8 / num_samples if num_samples > 0 else 0.25

    # Plot length distribution
    for i, (sample_name, df_unique_length) in enumerate(all_unique_data.items()):
        # Simplify sample names, remove suffixes
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

        # If the sample name contains "-KD", retain only that part.
        if '-KD' in simplified_name:
            parts = simplified_name.split('-KD')
            if len(parts) > 1:
                simplified_name = parts[0] + '-KD' + parts[1].split('_')[0]

        x = df_unique_length['Length'] + i * bar_width - (num_samples * bar_width / 2) + (bar_width / 2)
        plt.bar(x, df_unique_length['UniqueCount'], width=bar_width, color=colors[i % len(colors)], label=simplified_name)

    plt.xlabel('Length (nt)')
    plt.ylabel('Number of unique sequences')
    plt.title('Unique Reads Length Distribution Comparison')
    plt.legend()

    # Adjust layout
    plt.tight_layout()

    # Save figure, using simple file name, PDF format
    plt.savefig(f"{output_prefix}/unique_length_distribution.pdf", format='pdf')
    plt.close()

def extract_sample_name(file_path):
    """Extract simple sample name from file path"""
    # Get file name (without path)
    file_name = os.path.basename(file_path)

    # Remove file extension and common suffixes
    sample_name = file_name
    for suffix in ['.map', '.fa.rmdup.map', '.fa.collapsed.map', '.fa.collapsed.no-dust.map', '.collapsed.no-dust.map', '.no-dust.map', '.fa.collapsed.no-dust', '.collapsed.no-dust', '.no-dust', '.fa']:
        if sample_name.endswith(suffix):
            sample_name = sample_name[:-len(suffix)]
            break

    # Remove "_1.fq" and similar suffixes
    patterns_to_remove = ['_1.fq', '_2.fq', '.fq', '_1.fastq', '_2.fastq', '.fastq']
    for pattern in patterns_to_remove:
        if pattern in sample_name:
            sample_name = sample_name.replace(pattern, '')

    # If the sample name contains "-KD", only keep that part.
    if '-KD' in sample_name:
        parts = sample_name.split('-KD')
        if len(parts) > 1:
            sample_name = parts[0] + '-KD' + parts[1].split('_')[0]

    return sample_name

def parse_map_file_optimized(map_file, ref_pirna_file='/home/gyk/project/ld_pirna/data/bmo.v3.0.fa', chunk_size=10000):
    """
    Optimized version of parse_map_file function, using chunk processing and more efficient data structures

    Parameters:
        map_file: map file path
        ref_pirna_file: reference piRNA file path
        chunk_size: number of lines to process at a time

    Returns:
        read_info, read_to_pirna_count, read_to_pirnas, filter_stats
    """
    # Store the number of reference piRNAs each read aligns to
    read_to_pirna_count = defaultdict(int)
    # Store the count and length of each read
    read_info = {}
    # Store the reference piRNAs and positions each read aligns to
    read_to_pirnas = defaultdict(list)

    # Record the number of filtered sequences
    filtered_count = 0
    total_count = 0

    # Start timer
    start_time = time.time()

    print(f"Parsing reference piRNA file: {ref_pirna_file}")
    # Read sequences from reference piRNA file
    pirna_sequences = {}
    try:
        # Use iterator mode to reduce memory usage
        for record in SeqIO.parse(ref_pirna_file, "fasta"):
            pirna_sequences[record.id] = str(record.seq)
        print(f"Successfully read {len(pirna_sequences)} reference piRNA sequences, took {time.time() - start_time:.2f} seconds")
    except Exception as e:
        print(f"Warning: Error reading reference piRNA file: {e}")
        print("Will use target_seq as alternative")

    # Find the corresponding .fa.collapsed file to get sequence count information
    collapsed_file = map_file.replace('.map', '')
    if not os.path.exists(collapsed_file):
        # Try other possible file name patterns
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

    # Use Biopython to read sequence count information and sequences from collapsed file
    seq_counts = {}
    read_sequences = {}
    if os.path.exists(collapsed_file):
        collapsed_start_time = time.time()
        print(f"Parse the collapsed file:{collapsed_file}")
        try:
            #   Use tqdm to show processing progress
            with open(collapsed_file, 'r') as f:
                # Calculate file line count to set progress bar
                line_count = sum(1 for _ in f)

            record_count = 0
            for record in tqdm(SeqIO.parse(collapsed_file, "fasta"), total=line_count//2, desc="Read sequences"):
                # Extract count information from FASTA header
                count = 1
                if '-' in record.id:
                    try:
                        count = int(record.id.split('-')[-1])
                    except ValueError:
                        pass

                # Storage sequences and their counts
                seq = str(record.seq)
                seq_counts[seq] = count

                # Storage sequence ID and sequence mapping relationship
                read_sequences[record.id] = seq
                record_count += 1

            print(f"Successfully read {len(seq_counts)} sequences and their counts, took {time.time() - collapsed_start_time:.2f} seconds")
            print(f"Successfully stored {len(read_sequences)} sequence ID and sequence mapping relationships")
        except Exception as e:
            print(f"Warning: Error parsing collapsed file: {e}")
    else:
        print(f"Warning: collapsed file not found: {collapsed_file}")

    # Read map file
    print(f"Parsing map file: {map_file}")
    try:
        # Use pandas to read TSV format map file in chunks to reduce memory usage
        # First get column names
        header = pd.read_csv(map_file, sep='\t', nrows=1, header=None)
        column_count = len(header.columns)

        # Set column names
        column_names = ['trans_id', 'trans_coord', 'target_seq', 'probe_id', 'probe_seq', 'num_mismatch']
        if column_count > 6:
            column_names += [f'col{i+7}' for i in range(column_count-6)]

        # Calculate file line count to set progress bar
        with open(map_file, 'r') as f:
            line_count = sum(1 for _ in f) - 1  # Skip header row

        # Read in chunks to process
        reader = pd.read_csv(map_file, sep='\t', skiprows=1, header=None, names=column_names, chunksize=chunk_size)

        # Show progress bar
        with tqdm(total=line_count, desc="Parsing map file") as pbar:
            for chunk in reader:
                chunk_size = len(chunk)
                total_count += chunk_size
                pbar.update(chunk_size)

                for _, row in chunk.iterrows():
                    pirna_id = row['trans_id']  # Use trans_id as pirna_id
                    raw_start_pos = int(row['trans_coord'])

                    # Adjust start_pos, because using 3-23nt of each sequence for alignment
                    start_pos = raw_start_pos - 3

                    # Get ref_seq from reference piRNA file
                    if pirna_id in pirna_sequences:
                        ref_seq = pirna_sequences[pirna_id]
                    else:
                        # If not found, use target_seq as alternative
                        ref_seq = row['target_seq']
                        if total_count < 10:  # show only a few warnings
                            print(f"Warning: cannot find reference sequence for piRNA ID: {pirna_id}, using target_seq as alternative")

                    probe_id = row['probe_id']
                    read_seq = row['probe_seq']

                    # If probe_id in read_sequences, use mapped sequence first
                    if probe_id in read_sequences:
                        read_seq = read_sequences[probe_id]

                    _mismatches = int(row['num_mismatch'])

                    # Filter sequences with alignment start position > 6 (original position > 9)
                    if raw_start_pos > 9:  # Equivalent to adjusted start_pos > 6
                        filtered_count += 1
                        continue

                    # Extract count information from probe_id, format as ID-COUNT
                    read_count = 1  # default count = 1
                    if '-' in probe_id:
                        try:
                            read_count = int(probe_id.split('-')[-1])
                        except ValueError:
                            pass

                    # Get sequence count from seq_counts if available
                    if read_seq in seq_counts:
                        read_count = seq_counts[read_seq]

                    # Use read sequence as key, assume all are positive strand
                    strand = '+'
                    read_key = f"{read_seq}_{strand}"

                    # Update read sequence count and length
                    read_info[read_key] = (read_count, len(read_seq))

                    # Update read sequence alignment to reference piRNA
                    read_to_pirna_count[read_key] += 1
                    read_to_pirnas[read_key].append((pirna_id, start_pos, len(ref_seq)))
    except Exception as e:
        print(f"Warning: error parsing map file: {e}")
        import traceback
        traceback.print_exc()

    # Return filtering statistics for detailed mode display
    filter_stats = (total_count, filtered_count)
    print(f"Parsing complete: total {total_count} records, filtered {filtered_count}, elapsed {time.time() - start_time:.2f} s")
    return read_info, read_to_pirna_count, read_to_pirnas, filter_stats

def process_file(map_file, ref_file, output_prefix, verbose=False):
    """
    Function to process a single map file for parallel processing
    Only responsible for data processing and saving, not for plotting
    """
    sample_name = extract_sample_name(map_file)
    if verbose:
        print(f"Processing sample {sample_name}: {map_file}")

    # Use optimized parsing function
    read_info, read_to_pirna_count, read_to_pirnas, filter_stats = parse_map_file_optimized(map_file, ref_file)

    if verbose:
        total_reads, filtered_reads = filter_stats
        print(f"  Found {len(read_info)} unique reads in {sample_name}")
        print(f"  Filtered {filtered_reads}/{total_reads} reads with start position > 9 (adjusted > 6)")

    # Analyze length distribution
    df_length = analyze_length_distribution(read_info, read_to_pirna_count)

    # Analyze 5' end shift distribution
    df_5prime = analyze_5prime_shift(read_info, read_to_pirna_count, read_to_pirnas)

    # Analyze 3' end relative position distribution
    df_3prime = analyze_3prime_position(read_info, read_to_pirna_count, read_to_pirnas)

    # Analyze unique reads length distribution
    df_unique = analyze_unique_reads_length(read_info)

    # Create output directory
    os.makedirs(output_prefix, exist_ok=True)

    # Save results
    df_length.to_csv(f"{output_prefix}/{sample_name}_length_distribution.csv")
    df_5prime.to_csv(f"{output_prefix}/{sample_name}_5prime_shift.csv")
    df_3prime.to_csv(f"{output_prefix}/{sample_name}_3prime_position.csv")
    df_unique.to_csv(f"{output_prefix}/{sample_name}_unique_reads_length.csv")

    # Return sample name and processing results
    return sample_name

def main():
    # Configure CLI
    parser = argparse.ArgumentParser(description="Analyze piRNA length distribution, 5' end shift, 3' end relative position, and unique reads length distribution")
    parser.add_argument('-i', '--input', nargs='+', required=True, help='One or more input .map files')
    parser.add_argument('-o', '--output', required=True, help='Output directory or prefix')
    parser.add_argument('-v', '--verbose', action='store_true', help='Show verbose processing information')
    parser.add_argument('-p', '--prefix', default='', help='Output filename prefix (default: empty)')
    parser.add_argument('--no-biopython', action='store_true', help='Do not use Biopython; fall back to legacy parsing')
    parser.add_argument('-r', '--ref', default='/home/gyk/project/ld_pirna/data/bmo.v3.0.fa', help='Reference piRNA FASTA file path')
    parser.add_argument('--threads', type=int, default=0, help='Number of worker processes (0 uses all available CPU cores)')

    # Parse command line arguments
    args = parser.parse_args()

    # Get input file list and output prefix
    map_files = args.input
    output_prefix = args.output
    verbose = args.verbose
    file_prefix = args.prefix
    use_biopython = not args.no_biopython

    # Set number of worker processes
    threads = args.threads
    if threads <= 0:
        threads = mp.cpu_count()

    if use_biopython:
        try:
            import Bio
            if verbose:
                print(f"Using Biopython {Bio.__version__} for sequence analysis")
        except ImportError:
            print("Warning: Biopython not found; using legacy parsing")
            use_biopython = False

    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_prefix) if os.path.dirname(output_prefix) else output_prefix
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print(f"Using {threads} CPU core(s) for processing")

    # Process multiple samples
    all_data = {}
    all_unique_data = {}
    start_time = time.time()

    # Use multiprocessing to process files in parallel
    if len(map_files) > 1 and threads > 1:
        print(f"Processing {len(map_files)} files in parallel...")
        # Create process pool
        pool = mp.Pool(min(threads, len(map_files)))

        # Prepare process function arguments
        process_args = [(map_file, args.ref, output_prefix, verbose) for map_file in map_files]

        # Parallel process files
        results = []
        for result in tqdm(pool.starmap(process_file, process_args), total=len(map_files), desc="Processing files"):
            results.append(result)

        # Close process pool
        pool.close()
        pool.join()

        print("Data processing complete, start generating figures...")

        # Collect data from all samples for combined plots
        for map_file in map_files:
            sample_name = extract_sample_name(map_file)

            # Read saved data files
            df_length = pd.read_csv(f"{output_prefix}/{sample_name}_length_distribution.csv")
            df_shift = pd.read_csv(f"{output_prefix}/{sample_name}_5prime_shift.csv")
            df_position = pd.read_csv(f"{output_prefix}/{sample_name}_3prime_position.csv")
            df_unique_length = pd.read_csv(f"{output_prefix}/{sample_name}_unique_reads_length.csv")

            # Generate figures in main process
            if verbose:
                print(f"  Generating figures for {sample_name}...")

            # Plot various distributions
            try:
                # Plot with matplotlib
                plt.figure(figsize=(10, 6))
                plt.bar(df_length['Length'], df_length['Count'], color='#1f77b4')
                plt.xlabel('Length (nt)')
                plt.ylabel('Count')
                plt.title(f'{sample_name} - Length Distribution')
                plt.tight_layout()
                plt.savefig(f"{output_prefix}/{sample_name}_length_distribution.png")
                plt.close()

                # 5' end shift plot
                plt.figure(figsize=(10, 6))
                plt.bar(df_shift['Shift'], df_shift['Percentage'], color='#ff7f0e')
                plt.xlabel("5' position shift")
                plt.ylabel('% of analyzed piRNA loci')
                plt.title(f"{sample_name} - 5' End Shift")
                plt.tight_layout()
                plt.savefig(f"{output_prefix}/{sample_name}_5prime_shift.png")
                plt.close()

                # 3' end relative position plot
                plt.figure(figsize=(10, 6))
                plt.plot(df_position['Position'], df_position['Percentage'], 'o-', color='#2ca02c', linewidth=2)
                plt.xlabel("3' relative position")
                plt.ylabel('% of analyzed piRNA loci')
                plt.title(f"{sample_name} - 3' End Position")
                plt.tight_layout()
                plt.savefig(f"{output_prefix}/{sample_name}_3prime_position.png")
                plt.close()

                # Unique reads length distribution
                plt.figure(figsize=(10, 6))
                plt.bar(df_unique_length['Length'], df_unique_length['UniqueCount'], color='#d62728')
                plt.xlabel('Length (nt)')
                plt.ylabel('Number of unique sequences')
                plt.title(f'{sample_name} - Unique Reads Length Distribution')
                plt.tight_layout()
                plt.savefig(f"{output_prefix}/{sample_name}_unique_reads_length.png")
                plt.close()
            except Exception as e:
                print(f"Warning: error generating figures for {sample_name}: {e}")

            # Add prefix and rename files (if needed)
            if file_prefix:
                os.rename(f"{output_prefix}/{sample_name}_length_distribution.csv",
                          f"{output_prefix}/{file_prefix}{sample_name}_length_distribution.tsv")
                os.rename(f"{output_prefix}/{sample_name}_5prime_shift.csv",
                          f"{output_prefix}/{file_prefix}{sample_name}_5prime_shift.tsv")
                os.rename(f"{output_prefix}/{sample_name}_3prime_position.csv",
                          f"{output_prefix}/{file_prefix}{sample_name}_3prime_position.tsv")
                os.rename(f"{output_prefix}/{sample_name}_unique_reads_length.csv",
                          f"{output_prefix}/{file_prefix}{sample_name}_unique_length_distribution.tsv")

                # Rename image files
                os.rename(f"{output_prefix}/{sample_name}_length_distribution.png",
                          f"{output_prefix}/{file_prefix}{sample_name}_length_distribution.png")
                os.rename(f"{output_prefix}/{sample_name}_5prime_shift.png",
                          f"{output_prefix}/{file_prefix}{sample_name}_5prime_shift.png")
                os.rename(f"{output_prefix}/{sample_name}_3prime_position.png",
                          f"{output_prefix}/{file_prefix}{sample_name}_3prime_position.png")
                os.rename(f"{output_prefix}/{sample_name}_unique_reads_length.png",
                          f"{output_prefix}/{file_prefix}{sample_name}_unique_reads_length.png")

            # Store data for combined plots
            all_data[sample_name] = (df_length, df_shift, df_position)
            all_unique_data[sample_name] = df_unique_length
    else:
        # Single-process file handling
        for map_file in map_files:
            # Extract sample name from filename
            sample_name = extract_sample_name(map_file)

            if verbose:
                print(f"Processing sample {sample_name}: {map_file}")
            else:
                print(f"Processing sample {sample_name}...")

            # Use optimized parsing function
            read_info, read_to_pirna_count, read_to_pirnas, filter_stats = parse_map_file_optimized(map_file, args.ref)

            if verbose:
                total_reads, filtered_reads = filter_stats
                print(f"  Found {len(read_info)} unique reads in {sample_name}")
                print(f"  Filtered {filtered_reads}/{total_reads} reads with start position > 9 (adjusted > 6)")

            # Analyze length distribution
            df_length = analyze_length_distribution(read_info, read_to_pirna_count)

            # Analyze 5' end shift distribution
            df_shift = analyze_5prime_shift(read_info, read_to_pirna_count, read_to_pirnas)

            # Analyze 3' end relative position distribution
            df_position = analyze_3prime_position(read_info, read_to_pirna_count, read_to_pirnas)

            # Analyze unique reads length distribution
            df_unique_length = analyze_unique_reads_length(read_info)

            # Save data with concise filenames
            os.makedirs(f"{output_prefix}", exist_ok=True)
            df_length.to_csv(f"{output_prefix}/{file_prefix}{sample_name}_length_distribution.tsv", sep='\t', index=False)
            df_shift.to_csv(f"{output_prefix}/{file_prefix}{sample_name}_5prime_shift.tsv", sep='\t', index=False)
            df_position.to_csv(f"{output_prefix}/{file_prefix}{sample_name}_3prime_position.tsv", sep='\t', index=False)
            df_unique_length.to_csv(f"{output_prefix}/{file_prefix}{sample_name}_unique_length_distribution.tsv", sep='\t', index=False)

            if verbose:
                print(f"  Generating figures for {sample_name}...")

            # Drawing Charts
            try:
                # Using matplotlib to plot charts
                plt.figure(figsize=(10, 6))
                plt.bar(df_length['Length'], df_length['Count'], color='#1f77b4')
                plt.xlabel('Length (nt)')
                plt.ylabel('Count')
                plt.title(f'{sample_name} - Length Distribution')
                plt.tight_layout()
                plt.savefig(f"{output_prefix}/{file_prefix}{sample_name}_length_distribution.png")
                plt.close()

                # 5' end displacement distribution map
                plt.figure(figsize=(10, 6))
                plt.bar(df_shift['Shift'], df_shift['Percentage'], color='#ff7f0e')
                plt.xlabel("5' position shift")
                plt.ylabel('% of analyzed piRNA loci')
                plt.title(f"{sample_name} - 5' End Shift")
                plt.tight_layout()
                plt.savefig(f"{output_prefix}/{file_prefix}{sample_name}_5prime_shift.png")
                plt.close()

                # 3' end relative position distribution map
                plt.figure(figsize=(10, 6))
                plt.plot(df_position['Position'], df_position['Percentage'], 'o-', color='#2ca02c', linewidth=2)
                plt.xlabel("3' relative position")
                plt.ylabel('% of analyzed piRNA loci')
                plt.title(f"{sample_name} - 3' End Position")
                plt.tight_layout()
                plt.savefig(f"{output_prefix}/{file_prefix}{sample_name}_3prime_position.png")
                plt.close()

                # Unique reads length distribution map
                plt.figure(figsize=(10, 6))
                plt.bar(df_unique_length['Length'], df_unique_length['UniqueCount'], color='#d62728')
                plt.xlabel('Length (nt)')
                plt.ylabel('Number of unique sequences')
                plt.title(f'{sample_name} - Unique Reads Length Distribution')
                plt.tight_layout()
                plt.savefig(f"{output_prefix}/{file_prefix}{sample_name}_unique_reads_length.png")
                plt.close()
            except Exception as e:
                print(f"Error generating figures for {sample_name}: {e}")

            # Store data for combined plots
            all_data[sample_name] = (df_length, df_shift, df_position)
            all_unique_data[sample_name] = df_unique_length

    # Always plot combined distributions, even if only one sample
    if verbose:
        print("Generating combined distribution figures...")

    # Plot combined figures for regular distributions; export the three subplots as separate images
    plot_combined_distributions(all_data, f"{output_prefix}/{file_prefix}")

    # Plot combined figure for unique reads length distribution
    plot_combined_unique_reads_length(all_unique_data, f"{output_prefix}/{file_prefix}")

    print(f"Analysis complete. Results saved to: {output_prefix}")
    print(f"Total elapsed time: {time.time() - start_time:.2f} seconds")

if __name__ == "__main__":
    main()
