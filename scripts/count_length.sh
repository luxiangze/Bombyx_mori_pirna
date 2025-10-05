#!/bin/bash

# Input directory (contains *.fq.gz files)
target_dir=$1

# Output file
OUTPUT_FILE="$target_dir/length_distribution.txt"

# Create a temporary directory
TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

echo "Processing files in $target_dir..."

# Collect all files
FILES=("$target_dir"/*fq.gz)
if [ ${#FILES[@]} -eq 0 ]; then
    echo "No fq.gz files found in $target_dir"
    exit 1
fi

# Create or truncate the aggregated lengths file
> "$TMP_DIR/all_lengths.txt"

# For each file, create a temp file to store length counts
for file in "${FILES[@]}"; do
    base_filename=$(basename "$file" .fq.gz)
    echo "Processing $base_filename..."

    # Initialize an empty count file
    > "$TMP_DIR/$base_filename.counts"

    # Extract sequence length on every 2nd line and count occurrences
    zcat "$file" | awk 'NR % 4 == 2 {lengths[length($0)]++}
                         END {for (l in lengths) print l, lengths[l]}' > "$TMP_DIR/$base_filename.tmp"

    # Append lengths to global list
    cut -d' ' -f1 "$TMP_DIR/$base_filename.tmp" >> "$TMP_DIR/all_lengths.txt"

    # Collect counts
    while read -r len count; do
        echo "$len $count" >> "$TMP_DIR/$base_filename.counts"
    done < "$TMP_DIR/$base_filename.tmp"

    rm "$TMP_DIR/$base_filename.tmp"
done

# Get unique lengths and sort
sort -n -u "$TMP_DIR/all_lengths.txt" > "$TMP_DIR/unique_lengths.txt"

# Header line for output file
echo -n "Length" > "$OUTPUT_FILE"
for file in "${FILES[@]}"; do
    # Get base filename without suffix
    base_filename=$(basename "$file" .fq.gz)
    # Ensure .fq.gz suffix is removed
    clean_name=${base_filename%_fq.gz}
    echo -n -e "\t$clean_name" >> "$OUTPUT_FILE"
done
echo "" >> "$OUTPUT_FILE"

# Fill the table
echo "Creating the distribution table..."
while read -r length; do
    echo -n "$length" >> "$OUTPUT_FILE"

    for file in "${FILES[@]}"; do
        base_filename=$(basename "$file" .fq.gz)
        # Look up count for this length
        count=$(grep -w "^$length" "$TMP_DIR/$base_filename.counts" | awk '{print $2}' || echo 0)
        # Default to 0 if not found
        if [ -z "$count" ]; then
            count=0
        fi
        echo -n -e "\t$count" >> "$OUTPUT_FILE"
    done

    echo "" >> "$OUTPUT_FILE"
done < "$TMP_DIR/unique_lengths.txt"

echo "Done! Length distribution saved to $OUTPUT_FILE"