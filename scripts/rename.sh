#!/bin/bash

# Usage check
if [ $# -ne 2 ]; then
    echo "Usage: $0 <data_dir> <sample_map_file>"
    echo "Example: $0 data/raw_data data/sample_map.txt"
    exit 1
fi

# Args
DATA_DIR=$1
SAMPLE_MAP=$2

# Validate inputs
if [ ! -d "$DATA_DIR" ]; then
    echo "Error: directory '$DATA_DIR' does not exist"
    exit 1
fi

if [ ! -f "$SAMPLE_MAP" ]; then
    echo "Error: sample map file '$SAMPLE_MAP' does not exist"
    exit 1
fi

# Create md5 file
MD5_FILE="$DATA_DIR/md5.txt"
> "$MD5_FILE"

# Process mapping file
while read -r old_name new_name; do
  # Find matching files
  for file in "$DATA_DIR/${old_name}"*.fastq.gz; do
    if [[ -f "$file" ]]; then
      # Rename file
      new_file="$DATA_DIR/${new_name}.fq.gz"
      echo "Rename: $file -> $new_file"
      mv "$file" "$new_file"
      
      # Append new md5
      md5sum "$new_file" >> "$MD5_FILE"
    fi
  done
done < "$SAMPLE_MAP"

echo "All files renamed. MD5 checksums written to $MD5_FILE"