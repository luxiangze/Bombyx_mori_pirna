#!/bin/bash

# Script: remove_rRNA.sh
# Purpose: Remove rRNA sequences from small RNA data
# Usage: bash remove_rRNA.sh <input_dir> <bowtie_index_dir> <output_dir>

# Check arguments
if [ $# -ne 3 ]; then
    echo "Error: incorrect number of arguments"
    echo "Usage: bash $0 <input_dir> <bowtie_index_dir> <output_dir>"
    echo "Example: bash $0 /path/to/input /path/to/bowtie_index /path/to/output"
    exit 1
fi

# Get Parameters
INPUT_DIR="$1"
BOWTIE_INDEX_DIR="$2"
OUTPUT_DIR="$3"

# Create log directory
PROJECT_DIR=$(dirname $(dirname "$0"))
LOG_DIR="$PROJECT_DIR/logs"
mkdir -p "$LOG_DIR"

# Create log file (timestamped)
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="$LOG_DIR/remove_rRNA_${TIMESTAMP}.log"

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Start rRNA removal: $(date)" | tee -a "$LOG_FILE"
echo "Input dir: $INPUT_DIR" | tee -a "$LOG_FILE"
echo "Bowtie index dir: $BOWTIE_INDEX_DIR" | tee -a "$LOG_FILE"
echo "Output dir: $OUTPUT_DIR" | tee -a "$LOG_FILE"
echo "-----------------------------------" | tee -a "$LOG_FILE"

# Validate input dir
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: input dir '$INPUT_DIR' does not exist" | tee -a "$LOG_FILE"
    exit 1
fi

# Validate bowtie index dir
if [ ! -d "$BOWTIE_INDEX_DIR" ]; then
    echo "Error: bowtie index dir '$BOWTIE_INDEX_DIR' does not exist" | tee -a "$LOG_FILE"
    exit 1
fi

# Find all input files (support gz)
echo "Searching input files..." | tee -a "$LOG_FILE"
FILE_LIST=$(mktemp)
find "$INPUT_DIR" -type f \( -name "*.fa" -o -name "*.fasta" -o -name "*.fastq" -o -name "*.fq" -o -name "*.fa.gz" -o -name "*.fasta.gz" -o -name "*.fastq.gz" -o -name "*.fq.gz" \) > "$FILE_LIST"

# Count files
FILE_COUNT=$(wc -l < "$FILE_LIST")
echo "Found $FILE_COUNT files to process" | tee -a "$LOG_FILE"

# Detect CPU cores (use 90% of available)
TOTAL_CORES=$(nproc)
USABLE_CORES=$((TOTAL_CORES * 9 / 10))
if [ $USABLE_CORES -lt 1 ]; then
    USABLE_CORES=1
fi

# Compute per-file core allocation
if [ $FILE_COUNT -gt 0 ]; then
    # Max parallel jobs: min(file count, usable cores)
    MAX_PARALLEL_JOBS=$FILE_COUNT
    if [ $MAX_PARALLEL_JOBS -gt $USABLE_CORES ]; then
        MAX_PARALLEL_JOBS=$USABLE_CORES
    fi

    # Cores per file
    CORES_PER_FILE=1
    if [ $FILE_COUNT -lt $USABLE_CORES ]; then
        CORES_PER_FILE=$((USABLE_CORES / FILE_COUNT))
    fi

    echo "Parallel jobs: $MAX_PARALLEL_JOBS, cores per file: $CORES_PER_FILE" | tee -a "$LOG_FILE"
else
    echo "Warning: no files found" | tee -a "$LOG_FILE"
    rm -f "$FILE_LIST"
    exit 0
fi

# Define functions to handle individual files.
process_file() {
    local file=$1
    local cores=$2
    local log_file=$3

    # Determine filename and extension
    local filename=$(basename "$file")
    local is_gzipped=false
    local extension=""
    local basename=""

    # Check if gzipped
    if [[ "$filename" == *.gz ]]; then
        is_gzipped=true
        # Remove .gz extension to get original filename
        local orig_filename=${filename%.gz}
        extension="${orig_filename##*.}"
        basename="${orig_filename%.*}"
    else
        extension="${filename##*.}"
        basename="${filename%.*}"
    fi

    # Output path
    local output_file="$OUTPUT_DIR/${basename}.${extension}"

    # If input was gzipped, keep output gzipped
    if [ "$is_gzipped" = true ]; then
        output_file="${output_file}.gz"
    fi

    # Temp file for bowtie --un
    local temp_output="$(mktemp -p /tmp "${basename}_x_rRNA.XXXXXX")"

    # Log start
    echo "[PID $BASHPID] Start: $filename" >> "$log_file"

    # Choose input format option for bowtie
    local format_option=""
    if [[ "$extension" == "fa" || "$extension" == "fasta" ]]; then
        format_option="-f"
    elif [[ "$extension" == "fq" || "$extension" == "fastq" ]]; then
        format_option="-q"
    else
        echo "[PID $BASHPID]   Warning: unknown format, defaulting to FASTA" >> "$log_file"
        format_option="-f"
    fi

    # Bowtie can read gz directly; no need to decompress
    local input_for_bowtie="$file"

    # Run bowtie, capture unmapped (non-rRNA) reads
    echo "[PID $BASHPID]   Running bowtie..." >> "$log_file"
    bowtie -p $cores --un "$temp_output" $format_option \
           -v 2 -a --best --strata --norc \
           -S -x "$BOWTIE_INDEX_DIR/rRNA" \
           "$input_for_bowtie" \
           /dev/null 2>> "$log_file"

    # Check bowtie status
    if [ $? -eq 0 ]; then
        # Compress output if needed
        if [ "$is_gzipped" = true ]; then
            echo "[PID $BASHPID]   Compressing output..." >> "$log_file"
            gzip -c "$temp_output" > "$output_file"
        else
            mv "$temp_output" "$output_file"
        fi
        echo "[PID $BASHPID]   Success: non-rRNA saved to $output_file" >> "$log_file"
    else
        echo "[PID $BASHPID]   Error: failed processing $file" >> "$log_file"
    fi

    # Cleanup temp file
    if [ -f "$temp_output" ]; then
        rm -f "$temp_output"
    fi

    echo "[PID $BASHPID] Done: $filename" >> "$log_file"
    echo "[PID $BASHPID] -----------------------------------" >> "$log_file"
}

# Start processing
echo "Start processing..." | tee -a "$LOG_FILE"

# 初始化计数器
COMPLETED=0

# Read file list and dispatch jobs
while read -r file; do
    # Throttle concurrency to MAX_PARALLEL_JOBS
    while [ $(jobs -p | wc -l) -ge $MAX_PARALLEL_JOBS ]; do
        # 等待任一后台任务完成
        sleep 1
    done

    # Log dispatch
    filename=$(basename "$file")
    echo "Dispatch: $filename (progress: $COMPLETED/$FILE_COUNT)" | tee -a "$LOG_FILE"

    # Run in background
    process_file "$file" $CORES_PER_FILE "$LOG_FILE" &

    # Increase counter
    COMPLETED=$((COMPLETED + 1))
done < "$FILE_LIST"

# Wait for all jobs
echo "Waiting for all jobs to finish..." | tee -a "$LOG_FILE"
wait

# Cleanup file list
rm -f "$FILE_LIST"

echo "rRNA removal completed: $(date)" | tee -a "$LOG_FILE"
echo "Results saved to: $OUTPUT_DIR" | tee -a "$LOG_FILE"
echo "Log file: $LOG_FILE" | tee -a "$LOG_FILE"

exit 0
