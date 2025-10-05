#!/bin/bash

# Automated piRNA processing pipeline
# Usage: scripts/process_pirna.sh <input_dir> <output_dir> <reference_piRNA_fasta> [--skip-existing] [--reverse]

################################################################################
# Functions
################################################################################

# Logging helper
log() {
    local message="[$(date '+%Y-%m-%d %H:%M:%S')] $1"
    echo "$message" | tee -a "$MAIN_LOG"
}

# Simplify filename helper
simplify_filename() {
    local filename=$1
    # Remove common suffixes
    filename=$(echo "$filename" | sed -E 's/_1\.fq|_2\.fq|\.fq|_1\.fastq|_2\.fastq|\.fastq//g')
    echo "$filename"
}

# Process a single file
process_file() {
    local FASTQ_FILE=$1
    local CORES_TO_USE=6  # CPU cores allocated to this task
    local RAW_FILENAME=$(basename "$FASTQ_FILE" | sed 's/\.[^.]*$//')
    # Simplify filename
    local FILENAME=$(simplify_filename "$RAW_FILENAME")
    local FILE_LOG="$LOG_SUBDIR/${FILENAME}.log"

    # Skip if final output exists (check by simplified name to avoid conflicts)
    local FINAL_OUTPUT="$OUTPUT_DIR/${FILENAME}.fa.collapsed.map"
    if $SKIP_EXISTING && [ -f "$FINAL_OUTPUT" ]; then
        log "Skip processed file: $FILENAME (output exists)"
        return 0
    fi

    # Start
    log "Start processing: $FILENAME (PID: $$)"

    # 1) Convert fastq to fasta and pre-process
    log "[$FILENAME] Step 1: pre-processing sequences..." >> "$FILE_LOG" 2>&1
    log "[$FILENAME] Running seqkit with $CORES_TO_USE cores" >> "$FILE_LOG" 2>&1

    # Choose different processing workflows based on whether reverse processing is required.
    if $REVERSE_SEQ; then
        log "[$FILENAME] Reverse sequence mode enabled" >> "$FILE_LOG" 2>&1
        # Filter length, reverse sequences, then collapse
        seqkit fq2fa "$FASTQ_FILE" -j $CORES_TO_USE | \
        seqkit seq -m 23 -j $CORES_TO_USE | \
        seqkit seq -r -j $CORES_TO_USE | \
        fastx_collapser -o "$OUTPUT_DIR/${FILENAME}.fa.collapsed" >> "$FILE_LOG" 2>&1

        # Use pre-generated reversed reference
        local REF_DIR=$(dirname "$PIRNA_REF")
        local REVERSE_REF="$REF_DIR/$(basename "$PIRNA_REF").reverse"

        # 2) piRNA mapping (reversed reference)
        log "[$FILENAME] Step 2: mapping to reversed reference..." >> "$FILE_LOG" 2>&1
        # Use seqmap for mapping
        seqmap 0 "$OUTPUT_DIR/${FILENAME}.fa.collapsed" "$REVERSE_REF" "$OUTPUT_DIR/${FILENAME}.fa.collapsed.map" \
        /do_not_output_probe_without_match /forward_strand /cut:3,23 /output_all_matches >> "$FILE_LOG" 2>&1
        # Keep columns 4,2,1 (rename order) via awk
        awk '{print $4"\t"$2"\t"$1}' "$OUTPUT_DIR/${FILENAME}.fa.collapsed.map" > "$OUTPUT_DIR/${FILENAME}.fa.collapsed.map.filtered"
        mv "$OUTPUT_DIR/${FILENAME}.fa.collapsed.map.filtered" "$OUTPUT_DIR/${FILENAME}.fa.collapsed.map"
    else
        # Standard pipeline (no reverse)
        seqkit fq2fa "$FASTQ_FILE" -j $CORES_TO_USE | \
        seqkit seq -m 26 -j $CORES_TO_USE | \
        fastx_collapser -o "$OUTPUT_DIR/${FILENAME}.fa.collapsed" >> "$FILE_LOG" 2>&1

        # 2) piRNA mapping
        log "[$FILENAME] Step 2: mapping to reference..." >> "$FILE_LOG" 2>&1
        # Use seqmap for mapping
        seqmap 0 "$PIRNA_REF" "$OUTPUT_DIR/${FILENAME}.fa.collapsed" "$OUTPUT_DIR/${FILENAME}.fa.collapsed.map" \
        /do_not_output_probe_without_match /forward_strand /cut:3,23 /output_all_matches >> "$FILE_LOG" 2>&1
        # Reformat map: print 4, (3-2), 1; then filter where second column > -3
        awk '{print $4"\t"3-$2"\t"$1}' "$OUTPUT_DIR/${FILENAME}.fa.collapsed.map" | awk '$2 > -3' > "$OUTPUT_DIR/${FILENAME}.fa.collapsed.map.filtered"
        mv "$OUTPUT_DIR/${FILENAME}.fa.collapsed.map.filtered" "$OUTPUT_DIR/${FILENAME}.fa.collapsed.map"
    fi

    log "Done: $FILENAME"

    return 0
}

################################################################################
# Arguments and initialization
################################################################################

SCRIPT_PATH="$(cd "$(dirname "$0")" && pwd)" # absolute path to script dir
PROJECT_ROOT="$(dirname "$SCRIPT_PATH")" # project root is parent of scripts
LOG_DIR="$PROJECT_ROOT/logs"  # default logs dir
RESULTS_DIR="$PROJECT_ROOT/results"  # default results dir

# Check arguments
if [ $# -lt 4 ]; then
    echo "Usage: $0 <input_dir> <output_dir> <reference_piRNA_fasta> [--skip-existing] [--reverse]"
    exit 1
fi

# Obtain the necessary parameters
INPUT_DIR="$1"
OUTPUT_DIR="$2"
PIRNA_REF="$3"

shift 3  # shift positional args

# Parse flags --skip-existing / --reverse
SKIP_EXISTING=false
REVERSE_SEQ=false

# Process remaining arguments
while [ $# -gt 0 ]; do
    if [ "$1" = "--skip-existing" ]; then
        SKIP_EXISTING=true
    elif [ "$1" = "--reverse" ]; then
        REVERSE_SEQ=true
    fi
    shift
done

################################################################################
# Directories and files
################################################################################


TASK_ID=$(date +%Y%m%d_%H%M%S) # Get current time as task ID
mkdir -p "$OUTPUT_DIR" # Create output directory, log directory, and results directory
LOG_SUBDIR="$LOG_DIR/process_pirna_$TASK_ID"
mkdir -p "$LOG_SUBDIR"  # create dedicated log subdir
RESULTS_SUBDIR="$RESULTS_DIR/pirna_analysis_$TASK_ID"
mkdir -p "$RESULTS_SUBDIR" # create results dir
MAIN_LOG="$LOG_SUBDIR/process_pirna.log"
touch "$MAIN_LOG" # Create main log file


# Ensure pirna_analysis module is importable
export PYTHONPATH="$PROJECT_ROOT/scripts:$PYTHONPATH"

log "Start piRNA processing"
log "Script path: $SCRIPT_PATH"
log "Project root: $PROJECT_ROOT"
log "Input dir: $INPUT_DIR"
log "Output dir: $OUTPUT_DIR"
log "Log dir: $LOG_SUBDIR"
log "Results dir: $RESULTS_SUBDIR"
log "Reference: $PIRNA_REF"
log "Skip existing: $SKIP_EXISTING"
log "Reverse sequences: $REVERSE_SEQ"

# Check required tools
SCRIPTS_DIR="$SCRIPT_PATH"
TBR2_DIR="$SCRIPTS_DIR/TBr2"
if [ ! -d "$TBR2_DIR" ]; then
    log "Error: TBr2 scripts directory not found: $TBR2_DIR"
    exit 1
fi

# Check reference file
if [ ! -f "$PIRNA_REF" ]; then
    log "Error: reference file not found: $PIRNA_REF"
    exit 1
fi

# Check seqkit
if ! command -v seqkit &> /dev/null; then
    log "Error: seqkit not found in PATH"
    exit 1
fi

# Check seqmap
if ! command -v seqmap &> /dev/null; then
    log "Error: seqmap not found in PATH"
    exit 1
fi

# Check fastx_collapser
if ! command -v fastx_collapser &> /dev/null; then
    log "Error: fastx_collapser not found in PATH"
    exit 1
fi

################################################################################
# Reference analysis
################################################################################

# Basic stats on reference (once)
REF_STATS_FILE="$OUTPUT_DIR/$(basename "$PIRNA_REF").stats"
if $SKIP_EXISTING && [ -f "$REF_STATS_FILE" ]; then
    log "Skip reference stats: exists"
else
    log "Running reference basic analysis..."
    perl "$TBR2_DIR/TBr2_basic-analyses.pl" -i "$PIRNA_REF" -o "$REF_STATS_FILE" >> "$MAIN_LOG" 2>&1
fi

# Generate reversed reference if requested (once)
if $REVERSE_SEQ; then
    # Save reversed reference next to original
    REF_DIR=$(dirname "$PIRNA_REF")
    REVERSE_REF="$REF_DIR/$(basename "$PIRNA_REF").reverse"

    if [ -f "$REVERSE_REF" ]; then
        log "Skip reversed reference generation: exists ($REVERSE_REF)"
    else
        log "Generating reversed reference: $REVERSE_REF"
        log "Using 20 CPU cores for reverse generation"
        seqkit seq -r "$PIRNA_REF" -j 20 -o "$REVERSE_REF" >> "$MAIN_LOG" 2>&1

        if [ $? -eq 0 ]; then
            log "Reversed reference generated"
        else
            log "Error: failed to generate reversed reference"
            exit 1
        fi
    fi
fi

################################################################################
# File discovery and dispatch
################################################################################

# Discover fastq files
FASTQ_FILES=()
for ext in "fq.gz" "fastq.gz"; do
    if ls "$INPUT_DIR"/*.$ext 1> /dev/null 2>&1; then
        for file in "$INPUT_DIR"/*.$ext; do
            FASTQ_FILES+=("$file")
        done
    fi
done

# Ensure files found
if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
    log "Error: no fastq files found in input dir"
    exit 1
fi

log "Found ${#FASTQ_FILES[@]} files to process"

# Spawn background jobs
PIDS=()
for FASTQ_FILE in "${FASTQ_FILES[@]}"; do
    # Launch background process
    process_file "$FASTQ_FILE" &
    PIDS+=($!)
    log "Spawned PID $! for: $(basename "$FASTQ_FILE")"
done

# Wait for all processes
log "Waiting for all processes to finish..."
for pid in "${PIDS[@]}"; do
    wait $pid
    log "PID $pid finished"
done

log "All files processed!"

################################################################################
# Aggregate and analyze piRNA distributions
################################################################################

log "Start aggregating piRNA distributions..."
MAP_FILES=()
for FASTQ_FILE in "${FASTQ_FILES[@]}"; do
    # Original filename
    RAW_FILENAME=$(basename "$FASTQ_FILE" | sed 's/\.[^.]*$//')
    # Simplified filename
    SIMPLE_FILENAME=$(simplify_filename "$RAW_FILENAME")

    # Find map file by simplified filename
    MAP_FILE="$OUTPUT_DIR/${SIMPLE_FILENAME}.fa.collapsed.map"
    if [ -f "$MAP_FILE" ]; then
        MAP_FILES+=("$MAP_FILE")
    else
        log "Warning: map file not found: $MAP_FILE (raw: $RAW_FILENAME)"
    fi
done

if [ ${#MAP_FILES[@]} -gt 0 ]; then
    log "Found ${#MAP_FILES[@]} map files for aggregation"

    # Check whether the summary analysis results already exist.
    COMBINED_RESULT_EXISTS=false
    if $SKIP_EXISTING; then
        # Check if combined plots already exist
        COMBINED_PLOT="$RESULTS_SUBDIR/distributions.png"
        COMBINED_UNIQUE_PLOT="$RESULTS_SUBDIR/unique_length_distribution.png"

        if [ -f "$COMBINED_PLOT" ] && [ -f "$COMBINED_UNIQUE_PLOT" ]; then
            COMBINED_RESULT_EXISTS=true
            log "Skip aggregation: combined outputs exist"
        fi
    fi

    if ! $COMBINED_RESULT_EXISTS; then
        # Build CLI args list
        MAP_ARGS=""
        for MAP_FILE in "${MAP_FILES[@]}"; do
            # Append file paths after -i
            if [ -z "$MAP_ARGS" ]; then
                MAP_ARGS="-i $MAP_FILE"
            else
                MAP_ARGS="$MAP_ARGS $MAP_FILE"
            fi
        done

        # Threads equal to number of files (up to available cores)
        MAP_FILE_COUNT=${#MAP_FILES[@]}
        THREAD_ARGS="--threads $MAP_FILE_COUNT"
        ANALYSIS_CMD="python3 -m pirna_analysis.main $MAP_ARGS -o \"$RESULTS_SUBDIR\" -r $PIRNA_REF -v $THREAD_ARGS"
        log "Run: $ANALYSIS_CMD"

        # Execute aggregation
        eval $ANALYSIS_CMD >> "$MAIN_LOG" 2>&1
    fi

    if [ $? -eq 0 ]; then
        log "Aggregation completed"
    else
        log "Error: aggregation failed"
    fi
else
    log "Error: no map files found for aggregation"
fi

log "Detailed logs at: $LOG_SUBDIR"

# Write processing summary
SUMMARY_FILE="$LOG_SUBDIR/summary.txt"
{
    echo "========== piRNA Processing Summary ==========="
    echo "Time: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "Script: $SCRIPT_PATH"
    echo "Project root: $PROJECT_ROOT"
    echo "Input dir: $INPUT_DIR"
    echo "Output dir: $OUTPUT_DIR"
    echo "Log dir: $LOG_SUBDIR"
    echo "Results dir: $RESULTS_SUBDIR"
    echo "Reference: $PIRNA_REF"
    echo "Files processed: ${#FASTQ_FILES[@]}"
    echo "CPU cores per task: 10"
    echo "Skip existing: $SKIP_EXISTING"
    echo "Reverse mode: $REVERSE_SEQ"
    echo "==============================="
    echo ""
    echo "Files:"
    for file in "${FASTQ_FILES[@]}"; do
        echo "- $(basename "$file")"
    done
} > "$SUMMARY_FILE"

log "Summary saved to: $SUMMARY_FILE"

# Create simple README for results directory
README_FILE="$RESULTS_SUBDIR/README.txt"
{
    echo "========== piRNA Results README ==========="
    echo "Generated at: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "Input dir: $INPUT_DIR"
    echo "Files processed: ${#FASTQ_FILES[@]}"
    echo "Log dir: $LOG_SUBDIR"
    echo "Intermediate dir: $OUTPUT_DIR"
    echo "==============================="
    echo ""
    echo "Files:"
    echo "- combined_pirna_distribution_*.png: combined piRNA distribution plots"
    echo "- combined_pirna_distribution_*.tsv: combined piRNA distribution data"
    echo ""
    echo "Processed document:"
    for file in "${FASTQ_FILES[@]}"; do
        echo "- $(basename "$file")"
    done
} > "$README_FILE"

log "Results README saved to: $README_FILE"
