#!/bin/bash

# piPipes batch processing script
# Optimized for processing multiple fastq files through piPipes
input_dir="$1"
results_dir="$2"
log_dir="$3"
threads="${4:-70}"  # Default to 70 threads if not specified
genome="${5:-dm6}"  # Default genome is dm6 if not specified

# Validate input parameters
if [ -z "$input_dir" ] || [ -z "$results_dir" ] || [ -z "$log_dir" ]; then
    echo "Usage: $0 <input_dir> <results_dir> <log_dir> [threads] [genome]"
    echo "Example: $0 /path/to/fastq /path/to/results /path/to/logs 70 dm6"
    exit 1
fi

# Create output and log directories
mkdir -p "$results_dir"
mkdir -p "$log_dir"

echo "Starting piPipes batch processing..."
echo "Input directory: $input_dir"
echo "Results directory: $results_dir"
echo "Log directory: $log_dir"
echo "Using $threads threads and $genome genome"

# Count total files for progress reporting
total_files=$(ls "$input_dir"/*.fq.gz 2>/dev/null | wc -l)
current=0

for fq in "$input_dir"/*.fq.gz; do
    # Skip if no files found
    [ -e "$fq" ] || { echo "No .fq.gz files found in $input_dir"; exit 1; }
    
    sample=$(basename "$fq" .fq.gz)
    output_dir="$results_dir/${sample}.piPipes_out"
    current=$((current + 1))
    
    # Skip if output directory already exists
    if [ -d "$output_dir" ]; then
        echo "[$current/$total_files] Output directory $output_dir already exists, skipping..."
        continue
    fi
    
    echo "[$current/$total_files] Processing $sample..."
    
    # Run piPipes with logging
    piPipes small \
        -i "$fq" \
        -g "$genome" \
        -o "$output_dir" \
        -c "$threads" \
        -N uniqueXmiRNA \
        1> "$log_dir/${sample}.piPipes.stdout" \
        2> "$log_dir/${sample}.piPipes.stderr" &
        
    # Limit concurrent jobs to avoid system overload
    if [ $((current % 3)) -eq 0 ]; then
        wait
    fi
done

# Wait for all remaining jobs to complete
wait

echo "All processing complete. Results in $results_dir"