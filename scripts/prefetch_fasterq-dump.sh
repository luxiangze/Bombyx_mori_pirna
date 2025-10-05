#!/bin/bash

# Description: This script automatically downloads NCBI SRA reads
# Author: Asad Prodhan PhD
# Email: prodhan82@gmail.com
# Date: 2024-07-01
# Version: 1.0


# File containing the list of SRA accession numbers
SRA_LIST=$1
output_dir=$2

# Loop through each accession number in the list
while IFS= read -r accession; do
    echo "Processing $accession"
    prefetch $accession && fasterq-dump $accession --outdir reads
done < "$SRA_LIST"

# The end