#!/bin/bash

# Batch process fastq.gz/fq.gz files with fastp for QC and filtering
# Author: Yongkang Guo
# Created: 2025-03-24
# Updated: 2025-10-05

# Usage
usage() {
    echo "Usage: $0 [options] <input_dir> <output_dir>"
    echo "Options:"
    echo "  -p | --paired       Paired-end mode (default: single-end)"
    echo "  -h | --help         Show this help"
    exit 1
}

# Parse arguments
PAIRED_MODE=false

while [ $# -gt 0 ]; do
    case "$1" in
        -p|--paired)
            PAIRED_MODE=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            if [ -z "$INPUT_DIR" ]; then
                INPUT_DIR="$1"
            elif [ -z "$OUTPUT_DIR" ]; then
                OUTPUT_DIR="$1"
            else
                echo "Error: Unknown argument $1"
                usage
            fi
            shift
            ;;
    esac

# Validate required args
if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: <input_dir> and <output_dir> are required"
    usage
fi

REPORT_DIR="$OUTPUT_DIR/fastp_reports"

# Create output directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$REPORT_DIR"

# Info
if [ "$PAIRED_MODE" = true ]; then
    echo "Processing paired-end reads..."
else
    echo "Processing single-end reads (smRNA-seq mode)..."
fi

echo "Input dir: $INPUT_DIR"
echo "Output dir: $OUTPUT_DIR"
echo "Report dir: $REPORT_DIR"

# Counter
if [ "$PAIRED_MODE" = true ]; then
    # In paired-end mode, count number of R1 files
    total_files=$(find "$INPUT_DIR" \( -name "*_R1*.fastq.gz" -o -name "*_R1*.fq.gz" -o -name "*_1.fastq.gz" -o -name "*_1.fq.gz" \) | wc -l)
else
    # In single-end mode, count all files
    total_files=$(find "$INPUT_DIR" -name "*.fastq.gz" -o -name "*.fq.gz" | wc -l)
fi
processed=0

# Check files to process
if [ $total_files -eq 0 ]; then
    echo "No fastq.gz or fq.gz files found. Exiting."
    exit 1
fi

# Processing files according to schema
if [ "$PAIRED_MODE" = true ]; then
    # Paired-end: find all R1 files
    for r1_file in $(find "$INPUT_DIR" \( -name "*_R1*.fastq.gz" -o -name "*_R1*.fq.gz" -o -name "*_1.fastq.gz" -o -name "*_1.fq.gz" \)); do
        if [ ! -f "$r1_file" ]; then
            continue
        fi

        # Derive R2 filename pattern
        if [[ "$r1_file" == *_R1*.fastq.gz ]]; then
            r2_file="${r1_file/_R1/_R2}"
        elif [[ "$r1_file" == *_R1*.fq.gz ]]; then
            r2_file="${r1_file/_R1/_R2}"
        elif [[ "$r1_file" == *_1.fastq.gz ]]; then
            r2_file="${r1_file/_1/_2}"
        elif [[ "$r1_file" == *_1.fq.gz ]]; then
            r2_file="${r1_file/_1/_2}"
        fi

        # Ensure R2 exists
        if [ ! -f "$r2_file" ]; then
            echo "Warning: R2 not found: $r2_file, skipping $r1_file"
            continue
        fi

        # Base name
        if [[ "$r1_file" == *.fastq.gz ]]; then
            base_name=$(basename "$r1_file" | sed 's/_R1\(.*\)\.fastq\.gz$\|_1\(.*\)\.fastq\.gz$//')
        else
            base_name=$(basename "$r1_file" | sed 's/_R1\(.*\)\.fq\.gz$\|_1\(.*\)\.fq\.gz$//')
        fi

        echo "Processing paired sample: $base_name ($((processed+1)) / $total_files pairs)"

        # Output file paths
        output_r1="$OUTPUT_DIR/${base_name}_R1.clean.fastq.gz"
        output_r2="$OUTPUT_DIR/${base_name}_R2.clean.fastq.gz"

        # fastp reports
        html_report="$REPORT_DIR/${base_name}.html"
        json_report="$REPORT_DIR/${base_name}.json"

        # Run fastp (paired-end)
        fastp -i "$r1_file" -I "$r2_file" -o "$output_r1" -O "$output_r2" \
            --detect_adapter_for_pe --complexity_threshold 30 --thread 16 \
            --cut_mean_quality 20 --qualified_quality_phred 15 --unqualified_percent_limit 40 \
            --n_base_limit 5 --length_required 50 --low_complexity_filter \
            --html "$html_report" --json "$json_report"

        # Check fastp exit status
        if [ $? -eq 0 ]; then
            processed=$((processed + 1))
            echo "Done: $base_name"
        else
            echo "Failed: $base_name"
        fi

        echo "Progress: $processed / $total_files"
        echo "-----------------------------------"
    done
else
    # Single-end processing
    for fastq_file in "$INPUT_DIR"/*.fastq.gz "$INPUT_DIR"/*.fq.gz; do
        # Ensure file exists
        if [ ! -f "$fastq_file" ]; then
            continue
        fi

        # Get filename (without path and extension)
        if [[ "$fastq_file" == *.fastq.gz ]]; then
            filename=$(basename "$fastq_file" .fastq.gz)
        else
            filename=$(basename "$fastq_file" .fq.gz)
        fi

        # Output file path
        output_file="$OUTPUT_DIR/${filename}.clean.fastq.gz"

        # fastp report files
        html_report="$REPORT_DIR/${filename}.html"
        json_report="$REPORT_DIR/${filename}.json"

        echo "Processing: $filename"

        # Run fastp - smRNA-seq single-end parameters
        # Step 1: Remove adapter and initial filtering
        fastp -i "$fastq_file" -o "${OUTPUT_DIR}/${filename}.tmp1.fq.gz" \
            -a AACTGTAGGCACCATCAAT --length_required 1 --cut_front -W 4 -M 1 -q 1 -n 20 -w 16

        # Step 2: Additional filtering without adapter trimming
        fastp -i "${OUTPUT_DIR}/${filename}.tmp1.fq.gz" -o "${OUTPUT_DIR}/${filename}.tmp2.fq.gz" \
            -A -t 0 -w 16 --length_required 1 --cut_front -W 4 -M 1 -q 1 -n 20

        # Step 3: Final filtering and generate reports
        fastp -i "${OUTPUT_DIR}/${filename}.tmp2.fq.gz" -o "$output_file" \
            -A --length_required 18 --length_limit 37 -w 16 -n 0 -z 9 -q 20 -W 4 -M 20 \
            --html "$html_report" --json "$json_report"

        # Cleanup temp files
        rm "${OUTPUT_DIR}/${filename}.tmp1.fq.gz" "${OUTPUT_DIR}/${filename}.tmp2.fq.gz"


        # Check fastp exit status
        if [ $? -eq 0 ]; then
            processed=$((processed + 1))
            echo "Done: $filename"
        else
            echo "Failed: $filename"
        fi

        echo "Progress: $processed / $total_files"
        echo "-----------------------------------"
    done
fi

if [ "$PAIRED_MODE" = true ]; then
    echo "Completed. Processed $processed paired-end samples."
else
    echo "Completed. Processed $processed single-end samples."
fi
