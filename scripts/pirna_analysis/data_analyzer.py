#!/usr/bin/env python3
"""
Data analysis module - performs various analyses on parsed data.
"""

import pandas as pd
from collections import defaultdict

def analyze_length_distribution(read_info, read_to_pirna_count):
    """Analyze piRNA length distribution"""
    # Total normalized count per length
    length_counts = defaultdict(int)
    
    # Iterate over reads
    for read_key, (count, length) in read_info.items():
        # Number of reference piRNAs this read maps to
        pirna_count = read_to_pirna_count[read_key]
        # Normalized weight: read count / number of mapped reference piRNAs
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
    """Analyze 5' end shift distribution.
    
    In seqmap format, trans_coord is the alignment start position (1-based).
    """
    # Normalized counts per shift value
    shift_counts = defaultdict(float)
    
    # Iterate over reads
    for read_key, (count, read_length) in read_info.items():
        # Number of reference piRNAs this read maps to
        pirna_count = read_to_pirna_count[read_key]
        # Normalized weight: read count / number of mapped reference piRNAs
        normalized_weight = count / pirna_count
        
        # Update normalized counts per shift value
        for pirna_id, start_pos in read_to_pirnas[read_key]:
            shift = start_pos
            shift_counts[shift] += normalized_weight
    
    # Convert to DataFrame
    shifts = sorted(shift_counts.keys())
    counts = [shift_counts[s] for s in shifts]
    
    # Total
    total_count = sum(counts)
    
    # Percentage
    percentages = [count / total_count * 100 for count in counts] if total_count > 0 else [0] * len(counts)
    
    df_shift = pd.DataFrame({
        'Shift': shifts,
        'Count': counts,
        'Percentage': percentages
    })
    
    return df_shift

def analyze_3prime_position(read_info, read_to_pirna_count, read_to_pirnas):
    """Analyze 3' end relative position distribution"""
    # Normalized counts per 3' relative position
    position_counts = defaultdict(float)
    
    # Iterate over reads
    for read_key, (count, read_length) in read_info.items():
        # Number of reference piRNAs this read maps to
        pirna_count = read_to_pirna_count[read_key]
        # Normalized weight: read count / number of mapped reference piRNAs
        normalized_weight = count / pirna_count
        
        # Update counts per 3' relative position
        for pirna_id, start_pos in read_to_pirnas[read_key]:
            # Compute 3' end: start position + read length (seqmap start is 1-based)
            end_pos = start_pos + read_length
            position_counts[end_pos] += normalized_weight
    
    # Convert to DataFrame
    positions = sorted(position_counts.keys())
    counts = [position_counts[p] for p in positions]
    
    # Total
    total_count = sum(counts)
    
    # Percentage
    percentages = [count / total_count * 100 for count in counts] if total_count > 0 else [0] * len(counts)
    
    df_position = pd.DataFrame({
        'Position': positions,
        'Count': counts,
        'Percentage': percentages
    })
    
    return df_position

def analyze_unique_reads_length(read_info):
    """Analyze unique reads length distribution"""
    # Use pandas Series for efficient length counting
    lengths = [length for _, (_, length) in read_info.items()]
    length_series = pd.Series(lengths)
    
    # Count unique sequences per length
    unique_length_counts = length_series.value_counts().sort_index()
    
    # Convert to DataFrame
    df_unique_length = pd.DataFrame({
        'Length': unique_length_counts.index,
        'UniqueCount': unique_length_counts.values
    })
    
    return df_unique_length

def analyze_base_composition_by_position(read_info, read_to_pirna_count, read_to_pirnas, base_range=(-2, 46)):
    """
    Compute base composition percentage at each relative position.
    - base_range: (start, end), inclusive range of relative positions.
    Returns:
      DataFrame of base distribution (rows: relative positions; columns: A/T/C/G/N percentages)
    """
    from collections import Counter
    base_pos_dict = {pos: Counter() for pos in range(base_range[0], base_range[1]+1)}
    base_total = {pos: 0 for pos in range(base_range[0], base_range[1]+1)}
    base_types = set()

    # Accumulate base distribution per position
    for read_key, (count, read_length) in read_info.items():
        pirna_count = read_to_pirna_count[read_key]
        normalized_weight = count / pirna_count if pirna_count > 0 else 0
        seq = read_key
        for rel_pos in range(base_range[0], base_range[1]+1):
            seq_pos = rel_pos if rel_pos >= 0 else rel_pos + len(seq)
            if 0 <= seq_pos < len(seq):
                base = seq[seq_pos].upper()
                base_pos_dict[rel_pos][base] += normalized_weight
                base_total[rel_pos] += normalized_weight
                base_types.add(base)

    # Build percentage DataFrame
    base_types = sorted(list(base_types))
    base_dist_rows = []
    for rel_pos in range(base_range[0], base_range[1]+1):
        total = base_total[rel_pos]
        row = {'Position': rel_pos}
        for base in base_types:
            pct = (base_pos_dict[rel_pos][base] / total * 100) if total > 0 else 0
            row[base] = pct
        base_dist_rows.append(row)
    df_base_dist = pd.DataFrame(base_dist_rows)

    return df_base_dist

    
    
    
