#!/bin/bash

# Batch analyze piRNA trimming status
# Usage: ./batch_analyze_pirna.sh <map_dir> <output_dir> <sample_config>

# Check arguments
if [ $# -lt 3 ]; then
    echo "Usage: $0 <map_dir> <output_dir> <sample_config>"
    echo "Example: $0 /path/to/maps /path/to/output /path/to/config.txt"
    echo ""
    echo "Sample config format (space-separated 'sample type'):"
    echo "Control untreatment"
    echo "DHX15-KD3 treatment"
    echo "SUGP1-KD2 treatment"
    echo "RNPS1-KD1 treatment"
    exit 1
fi

# Args
MAP_DIR="$1"
OUTPUT_DIR="$2"
CONFIG_FILE="$3"
# Optional: log dir (default to output dir)
LOG_DIR="${4:-$OUTPUT_DIR}"

# Validate inputs
if [ ! -d "$MAP_DIR" ]; then
    echo "Error: directory $MAP_DIR does not exist"
    exit 1
fi

if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: config file $CONFIG_FILE does not exist"
    exit 1
fi

# Ensure output dir
mkdir -p "$OUTPUT_DIR"

# Read sample config; collect control/treatment
declare -A SAMPLE_TYPES
CONTROL_SAMPLES=()
TREATMENT_SAMPLES=()

while read -r sample type; do
    SAMPLE_TYPES["$sample"]="$type"
    if [ "$type" == "untreatment" ]; then
        CONTROL_SAMPLES+=("$sample")
    elif [ "$type" == "treatment" ]; then
        TREATMENT_SAMPLES+=("$sample")
    else
        echo "Warning: unknown sample type '$type', sample '$sample' will be ignored"
    fi
done < "$CONFIG_FILE"

# Ensure control/treatment found
if [ ${#CONTROL_SAMPLES[@]} -eq 0 ]; then
    echo "Error: no control (untreatment) samples found in config"
    exit 1
fi

if [ ${#TREATMENT_SAMPLES[@]} -eq 0 ]; then
    echo "Error: no treatment samples found in config"
    exit 1
fi

echo "Loaded from config:"
echo "  Control samples: ${CONTROL_SAMPLES[*]}"
echo "  Treatment samples: ${TREATMENT_SAMPLES[*]}"

# List all map files
declare -A MAP_FILES

echo "Map files found in $MAP_DIR:"
while IFS= read -r -d '' file; do
    echo "  $(basename "$file")"
done < <(find "$MAP_DIR" -name "*.map" -type f -print0)

# Match sample names to files
while IFS= read -r -d '' file; do
    filename=$(basename "$file")
    for sample in "${!SAMPLE_TYPES[@]}"; do
        if [[ "$filename" == *"$sample"* ]]; then
            MAP_FILES["$sample"]="$file"
            echo "Match: sample '$sample' -> file '$(basename "$file")'" 
            break
        fi
    done
done < <(find "$MAP_DIR" -name "*.map" -type f -print0)

# Check for missing sample files
MISSING_SAMPLES=()
for sample in "${!SAMPLE_TYPES[@]}"; do
    if [ -z "${MAP_FILES[$sample]}" ]; then
        MISSING_SAMPLES+=("$sample")
    fi
done

if [ ${#MISSING_SAMPLES[@]} -gt 0 ]; then
    echo "Warning: no map files found for the following samples:"
    for sample in "${MISSING_SAMPLES[@]}"; do
        echo "  $sample"
    done
    
    echo ""
    echo "Available map files:"
    while IFS= read -r -d '' file; do
        echo "  $(basename "$file")"
    done < <(find "$MAP_DIR" -name "*.map" -type f -print0)
    
    echo ""
    echo "Continue anyway? [y/N]"
    read -r answer
    if [[ ! "$answer" =~ ^[Yy]$ ]]; then
        echo "Analysis cancelled"
        exit 1
    fi
fi

# Timestamp and log setup
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
# Ensure log dir exists
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/analysis_log_$TIMESTAMP.log"

# For each control sample
for control_sample in "${CONTROL_SAMPLES[@]}"; do
    control_file="${MAP_FILES[$control_sample]}"
    
    # Skip if control file missing
    if [ -z "$control_file" ]; then
        echo "Warning: control sample $control_sample has no map file; skip"
        continue
    fi
    
    echo "Using $control_sample as control..." | tee -a "$LOG_FILE"
    
    # For each treatment sample
    for treatment_sample in "${TREATMENT_SAMPLES[@]}"; do
        treatment_file="${MAP_FILES[$treatment_sample]}"
        
        # Skip if treatment file missing
        if [ -z "$treatment_file" ]; then
            echo "  Warning: treatment sample $treatment_sample has no map file; skip" | tee -a "$LOG_FILE"
            continue
        fi
        
        echo "  Analyzing $control_sample vs $treatment_sample..." | tee -a "$LOG_FILE"
        
        # Create per-pair output dir
        sample_output_dir="$OUTPUT_DIR/${control_sample}_vs_${treatment_sample}"
        mkdir -p "$sample_output_dir"
        
        # Run analysis of 32nt and 28nt
        cmd="python $(dirname "$0")/identify_unprocessed_sequences.py -i \"$control_file\" \"$treatment_file\" -o \"$sample_output_dir\" -t both -v --control-name \"$control_sample\" --exp-name \"$treatment_sample\""
        echo "  $cmd" | tee -a "$LOG_FILE"
        
        # Execute
        python "$(dirname "$0")/identify_unprocessed_sequences.py" -i "$control_file" "$treatment_file" -o "$sample_output_dir" -t both -v --control-name "$control_sample" --exp-name "$treatment_sample" 2>&1 | tee -a "$LOG_FILE"
        
        # Normalize directory names if needed
        for dir in "$sample_output_dir"/*; do
            if [ -d "$dir" ]; then
                dir_name=$(basename "$dir")
                if [[ "$dir_name" == *"_vs_"* ]]; then
                    # Already concise name
                    continue
                fi
                # Create concise name
                new_dir_name="${control_sample}_vs_${treatment_sample}"
                if [ "$dir_name" != "$new_dir_name" ]; then
                    mv "$dir" "$sample_output_dir/$new_dir_name"
                fi
            fi
        done
        
        echo "  Completed: $control_sample vs $treatment_sample" | tee -a "$LOG_FILE"
        echo "" | tee -a "$LOG_FILE"
    done
done

echo "All analyses completed. Results at $OUTPUT_DIR"
echo "Log file: $LOG_FILE"
