#!/usr/bin/env Rscript

# Script to plot piRNA expression scatterplots
# Features:
# 1. Compute TPM for each piRNA
# 2. Scatterplot: x = control, y = treatment
# 3. Color points by p-value/FDR and logFC

# Check and install required R packages
required_packages <- c("ggplot2", "dplyr", "optparse", "viridis", "scales")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Load required packages
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(optparse)
  library(viridis)
  library(scales)
})

# Define CLI options
option_list <- list(
  make_option(c("-d", "--data_dir"), type = "character", default = NULL,
              help = "Directory containing edgeR result files"),
  make_option(c("-o", "--output_dir"), type = "character", default = NULL,
              help = "Output directory for figures"),
  make_option(c("-p", "--p_value_cutoff"), type = "numeric", default = 0.01,
              help = "Significance cutoff for P-value/FDR [default: %default]"),
  make_option(c("--use_fdr"), action = "store_true", default = FALSE,
              help = "Use FDR instead of P-value for significance"),
  make_option(c("--min_length"), type = "numeric", default = 0,
              help = "Minimum sequence length filter [default: %default]"),
  make_option(c("--max_length"), type = "numeric", default = 1000,
              help = "Maximum sequence length filter [default: %default]"),
  make_option(c("--pdf"), action = "store_true", default = FALSE,
              help = "Output PDF instead of PNG")
)

# Parse CLI options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$data_dir)) {
  stop("Please provide data directory path (-d or --data_dir)")
}

# If output_dir not provided, use current working directory
if (is.null(opt$output_dir)) {
  opt$output_dir <- getwd()
  message(paste("No output_dir provided; using current directory:", opt$output_dir))
}

# Ensure output directory exists
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# Function: compute TPM values
calculate_tpm <- function(counts) {
  # Ensure numeric vector
  counts <- as.numeric(counts)
  
  # Scaling factor per million reads
  scaling_factor <- sum(counts) / 1e6
  
  # Avoid division by zero
  if (scaling_factor == 0) {
    scaling_factor <- 1
  }
  
  # Compute TPM
  tpm <- counts / scaling_factor
  
  return(tpm)
}

# Process one file and create scatterplot
process_file <- function(file_path, output_dir, p_value_cutoff, output_pdf = FALSE, use_fdr = FALSE, min_length = 0, max_length = 1000) {
  # Read data
  message(paste("Processing file:", basename(file_path)))
  data <- read.csv(file_path, row.names = 1, check.names = FALSE)
  
  # Compute sequence lengths and filter
  seq_lengths <- nchar(rownames(data))
  message(paste("Sequences before filter:", nrow(data)))
  data_filtered <- data[seq_lengths >= min_length & seq_lengths <= max_length, ]
  message(paste("Sequences after filter:", nrow(data_filtered)))
  
  # If nothing remains after filtering, error out
  if (nrow(data_filtered) == 0) {
    stop("No sequences remain after length filtering")
  }
  
  # Use filtered data
  data <- data_filtered
  
  # Extract sample names
  file_name <- basename(file_path)
  sample_names <- strsplit(gsub("_edger_results_with_counts.csv", "", file_name), "_vs_")[[1]]
  control_name <- sample_names[1]
  treatment_name <- sample_names[2]
  
  # Compute TPM
  control_tpm <- calculate_tpm(data$control)
  treatment_tpm <- calculate_tpm(data$treatment)
  
  # Build plotting data
  plot_data <- data.frame(
    piRNA = rownames(data),
    control_tpm = control_tpm,
    treatment_tpm = treatment_tpm,
    logFC = data$logFC,
    PValue = data$PValue,
    FDR = data$FDR
  )
  
  # Choose P-value or FDR for significance
  if (use_fdr) {
    significance_value <- plot_data$FDR
    significance_type <- "FDR"
  } else {
    significance_value <- plot_data$PValue
    significance_type <- "P-value"
  }
  
  # Color categorization
  plot_data$color_category <- ifelse(significance_value <= p_value_cutoff, 
                                     ifelse(plot_data$logFC > 0, "Up", 
                                            ifelse(plot_data$logFC < 0, "Down", "Neutral")),
                                     "NonSig")
  
  # log10-transform TPM (add small constant to avoid log(0))
  plot_data$log10_control_tpm <- log10(plot_data$control_tpm + 1)
  plot_data$log10_treatment_tpm <- log10(plot_data$treatment_tpm + 1)
  
  # Limit logFC to [-6, 6]
  plot_data$logFC_limited <- pmax(pmin(plot_data$logFC, 6), -6)
  
  # Count up- and down-regulated points
  up_count <- sum(plot_data$color_category == "Up")
  down_count <- sum(plot_data$color_category == "Down")
  
  # Build scatter plot
  p <- ggplot(plot_data, aes(x = log10_control_tpm, y = log10_treatment_tpm)) +
    # Draw non-significant points first (gray)
    geom_point(data = subset(plot_data, color_category == "NonSig"), 
               color = "grey80", alpha = 0.4, size = 1) +
    # Then draw significant points (colored by logFC)
    geom_point(data = subset(plot_data, color_category != "NonSig"), 
               aes(color = logFC_limited), alpha = 0.4, size = 1.2) +
    # Add diagonal reference line
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
    # Color scale: blue-white-red
    scale_color_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      limits = c(-6, 6),
      name = "logFC",
      breaks = seq(-6, 6, by = 2),
      labels = seq(-6, 6, by = 2)
    ) +
    # Axis labels
    labs(
      x = paste("log10(TPM) -", control_name),
      y = paste("log10(TPM) -", treatment_name),
      title = paste(treatment_name, "vs", control_name, "piRNA Expression"),
      subtitle = paste(significance_type, "cutoff:", p_value_cutoff)
    ) +
    # Add annotations for counts
    annotate("text", x = min(plot_data$log10_control_tpm, na.rm = TRUE) + 1, 
             y = max(plot_data$log10_treatment_tpm, na.rm = TRUE) - 1, 
             label = paste("Up: n =", up_count), 
             color = "red", size = 4, hjust = 0) +
    annotate("text", x = max(plot_data$log10_control_tpm, na.rm = TRUE) - 1, 
             y = min(plot_data$log10_treatment_tpm, na.rm = TRUE) + 1, 
             label = paste("Down: n =", down_count), 
             color = "blue", size = 4, hjust = 1) +
    # Theme
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  # Save figure
  # Build filename (include length filter info)
  if (min_length > 0 || max_length < 1000) {
    length_suffix <- paste0("_len", min_length, "-", max_length)
  } else {
    length_suffix <- ""
  }
  
  output_file <- file.path(output_dir, 
                           paste0(treatment_name, "_vs_", control_name, 
                                  length_suffix,
                                  ifelse(output_pdf, ".pdf", ".png")))
  
  if (output_pdf) {
    ggsave(output_file, p, width = 8, height = 7, units = "in")
  } else {
    ggsave(output_file, p, width = 8, height = 7, units = "in", dpi = 300)
  }
  
  message(paste("Saved:", output_file))
  
  # 返回处理的数据
  return(list(plot = p, data = plot_data))
}

# Main
main <- function() {
  # List edgeR result files
  data_files <- list.files(opt$data_dir, 
                          pattern = "*_edger_results_with_counts.csv$", 
                          full.names = TRUE, 
                          recursive = FALSE)
  
  if (length(data_files) == 0) {
    stop(paste("No edgeR result files found in directory:", opt$data_dir))
  }
  
  message(paste("Found", length(data_files), "files to process"))
  
  # 处理每个文件
  results <- lapply(data_files, function(file) {
    tryCatch({
      process_file(file, opt$output_dir, opt$p_value_cutoff, opt$pdf, opt$use_fdr, opt$min_length, opt$max_length)
    }, error = function(e) {
      message(paste("Error processing file:", basename(file), "- Error:", e$message))
      NULL
    })
  })
  
  # 过滤掉NULL结果
  results <- results[!sapply(results, is.null)]
  
  if (length(results) > 0) {
    message(paste("Successfully processed", length(results), "files"))
  } else {
    message("No files were successfully processed")
  }
}

# 执行主函数
main()
