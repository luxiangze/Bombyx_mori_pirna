# piRNA Analysis Scripts Guide

This directory contains various scripts for analyzing piRNA data. These scripts are used to process sequencing data, analyze piRNA length distributions and positional features, and identify unprocessed sequences.

## Quick Start

Below is a complete example workflow for piRNA analysis:

```bash
# 1. Process raw sequencing data
./process_pirna.sh /path/to/fastq_files /path/to/output 8

# 2. Analyze piRNA length distribution and positional features
./batch_analyze_pirna.sh /path/to/output/map_files /path/to/results Control.map

# 3. Plot sample-to-sample comparison using the R script
Rscript plot_scatter_comparison.R --input_dir=/path/to/results --output_dir=/path/to/figures \
  --sample_config=/path/to/sample_config.txt
```

If you only need to replot figures from existing CSV results, use:

```bash
# Replot figures from previously generated CSV files
python3 replot_from_csv.py -i /path/to/csv_files -o /path/to/replotted_figures

# Or plot sample-to-sample comparison using the R script
Rscript plot_scatter_comparison.R --input_dir=/path/to/csv_files --output_dir=/path/to/figures
```

## Main Scripts

### Data Processing Scripts

#### `process_pirna.sh`

**Purpose**: Main driver to process piRNA data end-to-end, including preprocessing, mapping, and analysis.

**Usage**:
```bash
./process_pirna.sh <input_dir> <output_dir> [num_parallel] [reference_genome_path] [log_dir] [--skip-existing]
```

**Parameters**:
- `<input_dir>`: Directory containing raw FASTQ files
- `<output_dir>`: Output directory for processed results
- `[num_parallel]`: Optional, number of files to process in parallel (default: 4)
- `[reference_genome_path]`: Optional, path to piRNA reference genome file
- `[log_dir]`: Optional, directory to store logs
- `[--skip-existing]`: Optional, skip processing if results already exist

#### `process_fastq_with_fastp.sh`

**Purpose**: Process FASTQ files with fastp for quality control and filtering.

**Usage**:
```bash
./process_fastq_with_fastp.sh <input_path> <output_dir>
```

### piRNA Analysis Package

#### `pirna_analysis` package

This is a modular Python package for analyzing piRNA length distributions, 5' end shift distributions, and 3' end relative position distributions. It includes:

- `data_parser.py`: Parse map files and reference piRNA files
- `data_analyzer.py`: Analysis functions
- `visualization.py`: Plotting functions
- `utils.py`: Helper utilities
- `main.py`: CLI entry point and orchestration

#### `analyze_pirna.py`

**Purpose**: Entry-point script for piRNA analysis that calls the `pirna_analysis` package.

**Usage**:
```bash
python3 analyze_pirna.py -i <map_file1> [<map_file2> ...] -o <output_directory> [other options]
```

**Options**:
- `-i, --input`: One or more input .map files
- `-o, --output`: Output directory
- `-v, --verbose`: Verbose processing logs
- `-p, --prefix`: Output filename prefix
- `-r, --ref`: Reference piRNA file path
- `--threads`: Number of threads for parallel processing

#### `replot_from_csv.py`

**Purpose**: Replot piRNA analysis figures from previously generated CSV files to avoid re-parsing map files.

**Features**:
1. Auto-detect all CSV result files in the specified directory
2. Supports four plots: length distribution, 5' end shift distribution, 3' end relative position distribution, and unique read length distribution
3. Ensures integer x-axis tick marks for all plots
4. Supports custom output directory and file prefix

**Usage**:
```bash
python3 replot_from_csv.py -i <csv_dir> -o <output_dir> [-p <output_prefix>]
```

**Options**:
- `-i, --input`: Input directory containing CSV files
- `-o, --output`: Output directory
- `-p, --prefix`: Optional output file prefix

#### `analyze_pirna_distribution.py` (deprecated, modularized)

**Purpose**: Analyze piRNA length distribution, 5' end shift distribution, and 3' end relative position distribution.

**Note**: This script has been split into the `pirna_analysis` package. Use `analyze_pirna.py` or call `pirna_analysis.main` directly.

### Batch Scripts

#### `batch_analyze_pirna.sh`

**Purpose**: Batch-analyze all .map files in the specified directory with automatic CPU resource allocation.

**Features**:
1. Auto-discover all .map files under the given directory
2. Optional control file (default: Control.fa.collapsed.no-dust.map)
3. Dynamic CPU allocation targeting ~90% of server CPU usage
4. Analyze unprocessed sequences for both 32 nt and 28 nt
5. Timestamped log files
6. Detailed run information and error handling

### Sequence Analysis Scripts

#### `identify_unprocessed_sequences.py`

**Purpose**: Identify sequences not processed from 32 nt to 28 nt, and sequences not further processed from 28 nt to other lengths.

**Highlights**:
1. Focuses on a single control/treatment pair per run
2. Customizable sample names
3. Simplified output structure

**Usage**:
```bash
python3 identify_unprocessed_sequences.py <control.map> <treatment.map> -o <output_dir> [options]
```

#### `plot_length_distribution.py`

**Purpose**: Plot piRNA length distribution.

### Utility Scripts

#### `count_length.sh`

**Purpose**: Count sequence length distribution.

#### `rename.sh`

**Purpose**: Batch-rename files.

#### `piPipes_run.sh` and `tailor_run.sh`

**Purpose**: Invoke piPipes and Tailor tools for piRNA analysis.

## Examples

### Basic Processing Workflow

1. Process raw data with `process_pirna.sh`:
   ```bash
   ./process_pirna.sh /path/to/fastq_files /path/to/output 8
   ```

2. Analyze processed data with `analyze_pirna.py`:
   ```bash
   python3 analyze_pirna.py -i /path/to/output/*.map -o /path/to/results -v
   ```

3. Identify unprocessed sequences with `identify_unprocessed_sequences.py`:
   ```bash
   python3 identify_unprocessed_sequences.py /path/to/control.map /path/to/treatment.map -o /path/to/results
   ```

4. Replot figures from existing CSV files with `replot_from_csv.py`:
   ```bash
   python3 replot_from_csv.py -i /path/to/csv_files -o /path/to/replotted_figures -p "replotted_"
   ```

5. Plot sample comparisons with `plot_scatter_comparison.R`:
   ```bash
   Rscript plot_scatter_comparison.R --input_dir=/path/to/csv_files --output_dir=/path/to/output
   ```

### Batch Processing

Use `batch_analyze_pirna.sh` to analyze multiple samples in batch:
```bash
./batch_analyze_pirna.sh /path/to/map_files /path/to/output Control.map
```

### Data Visualization

#### `plot_scatter_comparison.R`

**Purpose**: General scatterplot tool to compare two samples, especially suitable for piRNA trimming data, using ggplot2 with extensive customization.

**Features**:
1. Auto-check and install required R packages
2. Handle multiple samples and generate pairwise comparisons
3. Dynamic color scale based on fold-change distribution
4. Support for log or linear axes
5. Clean visuals including diagonal reference line and colorbar
6. Support a sample config file to identify control/treatment
7. Auto-order axes with control on x-axis and treatment on y-axis

**Usage**:
```bash
Rscript plot_scatter_comparison.R --input_dir=/path/to/csv_files --output_dir=/path/to/output [--sample_config=/path/to/sample_config.txt]
```

**Options**:
- `--input_dir`: Directory of input CSV files (required)
- `--output_dir`: Directory to save output PDFs (required)
- `--sample_config`: Optional sample config file with sample names and types
- `--pattern`: File name pattern to select CSVs (default: "_pirna_trim.csv")
- `--id_column`: ID column for merging data (default: "piRNA_ID")
- `--value_column`: Numeric column used for comparison (default: "Trim_index")
- `--log_scale`: Whether to use log scales (default: TRUE)
- `--point_size`: Point size (default: 1)
- `--point_alpha`: Point transparency (default: 0.7)
- `--width`: PDF width in inches (default: 7)
- `--height`: PDF height in inches (default: 7)

**Sample config format**:
The sample config is a plain text file; each line contains:
```
<sample_name> <sample_type>
```
其中样本类型可以是`untreatment`（对照组）或`treatment`（处理组）。

Example:
```
Control untreatment
DHX15-KD3 treatment
SUGP1-KD2 treatment
RNPS1-KD1 treatment
```

**Examples**:
```bash
# Basic usage
Rscript plot_scatter_comparison.R --input_dir=/path/to/csv_files --output_dir=/path/to/output

# Use a sample config file
Rscript plot_scatter_comparison.R --input_dir=/path/to/csv_files --output_dir=/path/to/output \
  --sample_config=/path/to/sample_config.txt

# Custom parameters
Rscript plot_scatter_comparison.R --input_dir=/path/to/csv_files --output_dir=/path/to/output \
  --point_size=1.5 --point_alpha=0.5 --log_scale=FALSE --width=10 --height=8

# Use different file patterns and column names
Rscript plot_scatter_comparison.R --input_dir=/path/to/csv_files --output_dir=/path/to/output \
  --pattern="_expression.csv" --id_column="gene_id" --value_column="FPKM"
```

#### `plot_featrue_distribution.R`

**Purpose**: Draw donut plots to visualize category distributions using the tidyplots package.

**Features**:
1. Aesthetically pleasing donut plots via tidyplots
2. Customizable color palettes, styles, label sizes, etc.
3. Auto-calculate and display percentage per category
4. Support light and dark themes
5. Adjustable donut width and category order

**Usage**:
```bash
Rscript plot_featrue_distribution.R -i <CSV_path> -o <output_image_path> [options]
```

**Options**:
- `-i, --input`: Input CSV path, first column is category and second is value (required)
- `-o, --output`: Output image path (PDF/PNG etc.) (default: donut_plot.pdf)
- `-t, --title`: Plot title (default: "Feature distribution")
- `-w, --width`: Width in inches (default: 8)
- `-e, --height`: Height in inches (default: 8)
- `-c, --color_palette`: Palette name (friendly, neutral, vibrant, etc.) (default: friendly)
- `-l, --label_size`: Label text size (default: 12)
- `-p, --percentage`: Whether to show percentages (default: TRUE)
- `-d, --donut_width`: Donut width from 0 to 1 (default: 0.6)
- `-s, --style`: Plot style (light or dark) (default: light)
- `-r, --reverse`: Reverse category order (default: FALSE)

**Examples**:
```bash
# Basic usage
Rscript plot_featrue_distribution.R -i data/feature_counts.csv -o figures/feature_dist.pdf

# Custom parameters
Rscript plot_featrue_distribution.R -i data/feature_counts.csv -o figures/feature_dist.pdf \
  -t "Feature distribution" -c vibrant -s dark -d 0.8 -l 14 -r
```

**CSV format example**:
```
Category,Value
Promoter,150
Terminator,120
Intergenic,200
Exon,300
Intron,250
```

## Notes

1. Most scripts rely on specific file naming conventions; ensure inputs follow the expected patterns.
2. When processing large datasets, allocate CPU resources wisely to avoid overloading the server.
3. Use `--skip-existing` to skip previously processed files and save time.
4. R scripts (e.g., `plot_scatter_comparison.R`) require an R environment; the script auto-installs required packages.
5. For piRNA trimming comparisons, prefer the R script `plot_scatter_comparison.R` over the Python version to correctly render all points.
6. Large point counts (e.g., 400k+ per sample) can yield large PDFs (~8–10 MB); ensure sufficient disk space.

#### process_fastq_with_fastp.sh

```bash
# Single-end mode (for smRNA-seq)
scripts/process_fastq_with_fastp.sh <input_dir> <output_dir>

# Paired-end mode (for standard PE sequencing)
scripts/process_fastq_with_fastp.sh -p <input_dir> <output_dir>

# Help
scripts/process_fastq_with_fastp.sh -h