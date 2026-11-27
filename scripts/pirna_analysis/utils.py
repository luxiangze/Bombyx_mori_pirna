#!/usr/bin/env python3
"""
Utility functions module - provides common helpers.
"""

import os

def extract_sample_name(file_path):
    """Extract a concise sample name from a file path"""
    # Get basename without path
    file_name = os.path.basename(file_path)
    
    # Remove extension and common suffixes
    sample_name = file_name
    for suffix in ['.fa.rmdup.map', '.fa.collapsed.map', '.fa.collapsed.no-dust.map', '.collapsed.no-dust.map', '.no-dust.map', '.fa.collapsed.no-dust', '.collapsed.no-dust', '.no-dust', '.fa''.map']:
        if sample_name.endswith(suffix):
            sample_name = sample_name[:-len(suffix)]
            break
    
    # Remove sequencing-related suffix patterns
    patterns_to_remove = ['_1.fq', '_2.fq', '.fq', '_1.fastq', '_2.fastq', '.fastq']
    for pattern in patterns_to_remove:
        if pattern in sample_name:
            sample_name = sample_name.replace(pattern, '')
    
    # If name contains "-KD", keep only that segment
    if '-KD' in sample_name:
        parts = sample_name.split('-KD')
        if len(parts) > 1:
            sample_name = parts[0] + '-KD' + parts[1].split('_')[0]
            
    return sample_name
