#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert transposons predicted by HiTE into a GTF format compatible with TEtranscripts.
Usage: python hite_to_tetranscripts_gtf.py <input GFF> <output GTF>
"""

import sys
import os

def parse_attributes(attr_str):
    """Parse the attribute field in a GFF line"""
    attrs = {}
    for field in attr_str.split(';'):
        if not field.strip():
            continue
        key, value = field.split('=', 1)
        attrs[key] = value
    return attrs

def format_gtf_attributes(attrs):
    """Format attributes for the GTF attribute column"""
    result = []
    for key, value in attrs.items():
        result.append(f'{key} "{value}"')
    return '; '.join(result)

def convert_hite_to_tetranscripts(input_file, output_file):
    """Convert a HiTE GFF file to TEtranscripts-compatible GTF format"""
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            line = line.strip()
            # Skip comment lines
            if line.startswith('#'):
                continue

            # Parse a GFF line
            fields = line.split('\t')
            if len(fields) != 9:
                continue

            chrom, source, feature, start, end, score, strand, phase, attrs_str = fields

            # Parse attributes
            attrs = parse_attributes(attrs_str)

            # Read required attributes
            te_id = attrs.get('id', '')
            te_name = attrs.get('name', '')
            classification = attrs.get('classification', 'Unknown')

            # Handle classification info
            if '/' in classification:
                class_id, family_id = classification.split('/', 1)
            else:
                # Handle Unknown classification
                class_id = classification
                family_id = "Unknown"

            # Build new GTF attributes
            new_attrs = {
                "gene_id": te_name,
                "transcript_id": f"{te_id}",
                "class_id": class_id,
                "family_id": family_id
            }

            # Compose the GTF line
            new_line = [
                chrom,
                "HiTE",  # Set source
                "exon",  # Force all features to exon
                start,
                end,
                score,
                strand,
                phase,
                format_gtf_attributes(new_attrs)
            ]

            outfile.write('\t'.join(new_line) + '\n')

    print(f"Conversion completed. Output file: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input GFF> <output GTF>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    if not os.path.exists(input_file):
        print(f"Error: input file '{input_file}' does not exist!")
        sys.exit(1)

    convert_hite_to_tetranscripts(input_file, output_file)
