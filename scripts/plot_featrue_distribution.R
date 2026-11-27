#!/usr/bin/env Rscript

# Check and install required packages
packages <- c("tidyverse", "optparse", "tidyplots")
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)

# Load required packages
library(tidyverse)
library(optparse)
library(tidyplots)

# Parse CLI options
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input CSV path: first column=category, second column=value", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="donut_plot.pdf",
              help="Output image path (PDF/PNG) [default= %default]", metavar="character")
)

# Create option parser
opt_parser <- OptionParser(option_list=option_list,
                          description="\nDraw donut plots (tidyplots) from CSV.\nInput format: first column=category name, second column=value.")

# Parse options
opt <- parse_args(opt_parser)

# Validate input file argument
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input file is required. Use -i or --input.", call.=FALSE)
}

# Ensure input file exists
if (!file.exists(opt$input)) {
  stop(paste("Input file does not exist:", opt$input), call.=FALSE)
}

# Read data
data <- read_csv(opt$input, col_names = TRUE)

# Transform data for tidyplots
data_tidy <- data %>%
  pivot_longer(
    cols = -key1, # all columns except key1
    names_to = "Group", # move column names into new variable
    values_to = "Expression" # move values into variable
  )

# Create donut plot with tidyplots
data_tidy |>
  # base tidyplot with color mapping
  tidyplot(y = Expression, color = key1) |>
  # add donut layer
  add_donut() |>
  adjust_size(width = 250, height = 250) |>
  split_plot(by = Group) |>
  # save image
  save_plot(opt$output)

# Output info
cat(paste("\nDonut plot saved to:", normalizePath(opt$output), "\n"))
cat(paste("Number of categories:", nrow(data), "\n"))
