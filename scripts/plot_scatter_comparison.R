#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
#
# General scatterplot comparison tool
# Plot pairwise comparisons between two samples with customizable options
#
# Usage:
# Rscript plot_scatter_comparison.R [options]
#
# Options:
#   --input_dir=<dir>        Directory containing input CSV files
#   --output_dir=<dir>       Directory to save output PDF files
#   --pattern=<pattern>      Pattern for selecting CSV files [default: "_pirna_trim.csv"]
#   --id_column=<name>       ID column for merging [default: "piRNA_ID"]
#   --value_column=<name>    Numeric column to compare [default: "Trim_index"]
#   --log_scale=<TRUE/FALSE> Whether to use log axes [default: TRUE]
#   --point_size=<num>       Point size [default: 1]
#   --point_alpha=<num>      Point transparency [default: 0.7]
#   --width=<num>            PDF width (inches) [default: 7]
#   --height=<num>           PDF height (inches) [default: 7]
#   --help                   Show help

# Check and install missing packages
packages <- c("ggplot2", "dplyr", "tidyr", "readr", "optparse", "scales", "viridis")
missing_packages <- packages[!packages %in% installed.packages()[,"Package"]]
if(length(missing_packages) > 0) {
  message(paste0("Installing missing packages: ", paste(missing_packages, collapse=", ")))
  install.packages(missing_packages, repos="https://cloud.r-project.org")
}

# Load required packages
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(optparse)
  library(scales)
  library(viridis)
})

# Define CLI options
option_list <- list(
  make_option("--input_dir", type="character", default=NULL,
              help="Directory containing input CSV files"),
  make_option("--output_dir", type="character", default=NULL,
              help="Directory to save output PDF files"),
  make_option("--pattern", type="character", default="_pirna_trim.csv",
              help="Filename pattern to select CSV files [default: %default]"),
  make_option("--id_column", type="character", default="piRNA_ID",
              help="ID column to merge on [default: %default]"),
  make_option("--value_column", type="character", default="Trim_index",
              help="Numeric column to compare [default: %default]"),
  make_option("--log_scale", type="logical", default=TRUE,
              help="Use log axes [default: %default]"),
  make_option("--point_size", type="double", default=1,
              help="Point size [default: %default]"),
  make_option("--point_alpha", type="double", default=0.7,
              help="Point transparency [default: %default]"),
  make_option("--width", type="double", default=7,
              help="Output PDF width (inches) [default: %default]"),
  make_option("--height", type="double", default=7,
              help="Output PDF height (inches) [default: %default]"),
  make_option("--sample_config", type="character", default=NULL,
              help="Path to sample config (sample name and type)")
)

# Parse CLI options
opt_parser <- OptionParser(option_list=option_list,
                          description="General scatterplot comparison tool for two samples")
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$input_dir) || is.null(opt$output_dir)) {
  print_help(opt_parser)
  stop("Must specify both input_dir and output_dir", call.=FALSE)
}

# Create output directory
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# Simplify sample names
simplify_sample_name <- function(sample_name) {
  # Remove common suffixes
  suffixes <- c(".map", ".fa.rmdup.map", ".fa.collapsed.map", ".fa.collapsed.no-dust.map",
                ".collapsed.no-dust.map", ".no-dust.map", ".fa.collapsed.no-dust",
                ".collapsed.no-dust", ".no-dust", ".fa")

  for (suffix in suffixes) {
    if (endsWith(sample_name, suffix)) {
      sample_name <- substr(sample_name, 1, nchar(sample_name) - nchar(suffix))
      break
    }
  }

  # Remove sequencing-related suffixes
  patterns <- c("_1.fq", "_2.fq", ".fq", "_1.fastq", "_2.fastq", ".fastq")
  for (pattern in patterns) {
    sample_name <- gsub(pattern, "", sample_name)
  }

  # If sample name contains "-KD", keep only that part
  if (grepl("-KD", sample_name)) {
    parts <- strsplit(sample_name, "-KD")[[1]]
    if (length(parts) > 1) {
      kd_parts <- strsplit(parts[2], "_")[[1]]
      sample_name <- paste0(parts[1], "-KD", kd_parts[1])
    }
  }

  return(sample_name)
}

# List input files
input_files <- list.files(path = opt$input_dir, pattern = opt$pattern, full.names = TRUE)
if (length(input_files) < 2) {
  stop("At least two samples are required for comparison", call.=FALSE)
}

# Read sample config
sample_info <- list()
if (!is.null(opt$sample_config) && file.exists(opt$sample_config)) {
  message(paste("Reading sample config:", opt$sample_config))
  config_lines <- readLines(opt$sample_config)
  for (line in config_lines) {
    if (line != "") {
      parts <- strsplit(line, " ")[[1]]
      if (length(parts) >= 2) {
        sample_name <- parts[1]
        sample_type <- parts[2]
        sample_info[[sample_name]] <- sample_type
        message(paste("  Sample:", sample_name, ", Type:", sample_type))
      }
    }
  }
  message(paste("Sample config contains", length(sample_info), "samples"))
}

# Read all sample data
all_data <- list()
for (file in input_files) {
  # Extract sample name from filename
  file_name <- basename(file)
  sample_name <- gsub(opt$pattern, "", file_name)

  # Simplify sample name
  simple_name <- simplify_sample_name(sample_name)

  # Read CSV
  message(paste("Reading file:", file_name))
  df <- read_csv(file, show_col_types = FALSE)

  # Check mandatory columns
  if (!opt$id_column %in% colnames(df) || !opt$value_column %in% colnames(df)) {
    warning(paste("File", file_name, "missing required columns; skipped"))
    next
  }

  # Store data
  all_data[[simple_name]] <- df
  message(paste("  Rows:", nrow(df)))
  message(paste("  Columns:", paste(colnames(df), collapse=", ")))
}

message(paste("\nTotal samples:", length(all_data)))

# Plot pairwise comparisons
message("\nStart plotting sample comparisons...")
sample_names <- names(all_data)
n_samples <- length(sample_names)

# Helper: get sample type
get_sample_type <- function(sample_name) {
  # Try exact match first
  if (sample_name %in% names(sample_info)) {
    return(sample_info[[sample_name]])
  }

  # Else try partial match
  for (name in names(sample_info)) {
    if (grepl(name, sample_name, fixed = TRUE) || grepl(sample_name, name, fixed = TRUE)) {
      return(sample_info[[name]])
    }
  }

  # Fallback
  return("unknown")
}

# Decide sample order: control on x-axis, treatment on y-axis
should_swap_samples <- function(sample1, sample2) {
  type1 <- get_sample_type(sample1)
  type2 <- get_sample_type(sample2)

  # Swap if first is treatment and second is control
  if (type1 == "treatment" && type2 == "untreatment") {
    return(TRUE)
  }

  # Otherwise keep order
  return(FALSE)
}

# Check whether the pair is control vs treatment
is_control_vs_treatment <- function(sample1, sample2) {
  type1 <- get_sample_type(sample1)
  type2 <- get_sample_type(sample2)

  # TRUE when one is control and the other is treatment
  if ((type1 == "untreatment" && type2 == "treatment") ||
      (type1 == "treatment" && type2 == "untreatment")) {
    return(TRUE)
  }

  # Otherwise FALSE
  return(FALSE)
}

for (i in 1:(n_samples-1)) {
  for (j in (i+1):n_samples) {
    sample1_name <- sample_names[i]
    sample2_name <- sample_names[j]

    # Ensure it's a control vs treatment comparison
    if (!is_control_vs_treatment(sample1_name, sample2_name)) {
      message(paste("  Skipped:", sample1_name, "vs", sample2_name, "(not control vs treatment)"))
      next
    }

    # Swap order if needed
    if (should_swap_samples(sample1_name, sample2_name)) {
      temp <- sample1_name
      sample1_name <- sample2_name
      sample2_name <- temp
    }

    message(paste("  Plot:", sample1_name, "vs", sample2_name))

    # Sample types
    type1 <- get_sample_type(sample1_name)
    type2 <- get_sample_type(sample2_name)

    # Fetch data
    df1 <- all_data[[sample1_name]]
    df2 <- all_data[[sample2_name]]

    # Merge data
    merged_df <- df1 %>%
      select(!!sym(opt$id_column), !!sym(opt$value_column)) %>%
      inner_join(
        df2 %>% select(!!sym(opt$id_column), !!sym(opt$value_column)),
        by = opt$id_column,
        suffix = c("_1", "_2")
      )

    # Compute fold change
    merged_df <- merged_df %>%
      mutate(fold_change = .data[[paste0(opt$value_column, "_2")]] / (.data[[paste0(opt$value_column, "_1")]] + 1e-10))

    # Determine fold change range for color mapping
    fc_min <- min(merged_df$fold_change)
    fc_max <- max(merged_df$fold_change)

    # Ensure symmetry around 1
    if (fc_min < 1 && 1/fc_min > fc_max) {
      vmax <- 1/fc_min
      vmin <- fc_min
    } else {
      vmax <- fc_max
      vmin <- 1/fc_max
    }

    # Log fold change range for reference
    message(paste("    Fold change range:", round(fc_min, 2), "-", round(fc_max, 2)))

    # Color scale: dark blue - yellow - dark red
    color_scale <- scale_color_gradient2(
      low = "darkblue",
      mid = "yellow",
      high = "darkred",
      midpoint = 1,
      trans = "log10",
      name = "Fold change"
    )

    # Data range
    min_val <- min(c(merged_df[[paste0(opt$value_column, "_1")]], merged_df[[paste0(opt$value_column, "_2")]]))
    max_val <- max(c(merged_df[[paste0(opt$value_column, "_1")]], merged_df[[paste0(opt$value_column, "_2")]]))

    # Ensure positive min for log axes
    min_val <- max(min_val, 1e-10)

    # Build scatter plot
    p <- ggplot(merged_df, aes(
      x = .data[[paste0(opt$value_column, "_1")]],
      y = .data[[paste0(opt$value_column, "_2")]],
      color = fold_change
    )) +
      geom_point(size = opt$point_size, alpha = 0.4) + # point alpha = 0.4
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
      color_scale +
      labs(
        x = paste0(sample1_name, "\n(Normalized reads)"),
        y = paste0(sample2_name, "\n(Normalized reads)"),
        title = paste(sample1_name, "vs", sample2_name)
      ) +
      theme_bw() +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)
      )

    # Apply log axes if requested
    if (opt$log_scale) {
      p <- p +
        scale_x_log10(
          limits = c(min_val, max_val * 1.1),
          labels = function(x) format(x, scientific = FALSE, digits = 1)
        ) +
        scale_y_log10(
          limits = c(min_val, max_val * 1.1),
          labels = function(x) format(x, scientific = FALSE, digits = 1) # Using the correct logarithmic scale format.
        )
    } else {
      p <- p +
        scale_x_continuous(limits = c(0, max_val * 1.1)) +
        scale_y_continuous(limits = c(0, max_val * 1.1))
    }

    # No sample-type annotations

    # Save figure
    output_file <- file.path(opt$output_dir, paste0(sample1_name, "_vs_", sample2_name, ".pdf"))
    ggsave(output_file, plot = p, width = opt$width, height = opt$height)
    message(paste("    Saved:", output_file))
  }
}

message(paste("\nDone. Output saved to:", opt$output_dir))
