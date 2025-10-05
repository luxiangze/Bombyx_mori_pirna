#!/usr/bin/env Rscript

# TE differential expression analysis script (edgeR)
# Usage: Rscript te_differential_analysis.R [counts_matrix] [sample_config] [output_dir]

# Load required libraries
suppressPackageStartupMessages({
  if (!require("edgeR")) {
    install.packages("BiocManager", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
    BiocManager::install("edgeR")
    library(edgeR)
  }
  if (!require("ggplot2")) {
    install.packages("ggplot2", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
    library(ggplot2)
  }
})

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  counts_file <- "/home/gyk/project/ld_pirna/work/trasposon_analysis/TEcount/merged_TEcount.tsv"
  sample_config <- "/home/gyk/project/ld_pirna/data/sample_config.txt"
  output_dir <- "/home/gyk/project/ld_pirna/work/trasposon_analysis/TEcount/edgeR_results"
} else {
  counts_file <- args[1]
  sample_config <- args[2]
  output_dir <- args[3]
}

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Read sample configuration
sample_info <- read.table(sample_config, header = FALSE, stringsAsFactors = FALSE)
colnames(sample_info) <- c("sample_name", "condition")

# Read expression count data
counts_data <- read.delim(counts_file, header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

# Adjust sample name mapping for current data
# Note: mapping is hard-coded based on current column names
# If there are replicates, adjust the mapping accordingly
sample_mapping <- c(
  "control" = "Control",
  "DHX15" = "DHX15-KD3",
  "RNPS1" = "RNPS1-KD1",
  "SUGP1-Rep1" = "SUGP1-KD2",
  "SUGP1-Rep2" = "SUGP1-KD2"  # Note: Rep2 is mapped to the same sample name
)

# Align counts_data column names with sample_config
colnames_original <- colnames(counts_data)
colnames_mapped <- sapply(colnames_original, function(x) {
  if (x %in% names(sample_mapping)) {
    return(sample_mapping[x])
  } else {
    return(x)
  }
})
colnames(counts_data) <- colnames_mapped

# Subset samples present in both config and counts
sample_subset <- subset(sample_info, sample_name %in% colnames_mapped)

# Handle replicates (e.g., average multiple replicates)
if (length(unique(colnames_mapped)) < length(colnames_mapped)) {
  cat("Detected replicated samples, merging...\n")
  counts_merged <- matrix(0, nrow = nrow(counts_data), ncol = length(unique(colnames_mapped)))
  rownames(counts_merged) <- rownames(counts_data)
  colnames(counts_merged) <- unique(colnames_mapped)

  # 对重复样本进行合并（取平均值）
  for (col_name in unique(colnames_mapped)) {
    col_indices <- which(colnames_mapped == col_name)
    if (length(col_indices) > 1) {
      counts_merged[, col_name] <- rowMeans(counts_data[, col_indices, drop = FALSE])
    } else {
      counts_merged[, col_name] <- counts_data[, col_indices]
    }
  }
  counts_data <- counts_merged
}

# Ensure column order matches sample_subset
common_samples <- intersect(sample_subset$sample_name, colnames(counts_data))
counts_data <- counts_data[, common_samples]
sample_subset <- sample_subset[match(common_samples, sample_subset$sample_name), ]

# Create DGEList object
group <- factor(sample_subset$condition)
y <- DGEList(counts = counts_data, group = group)

# Filter lowly expressed TEs
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]

# Normalize
y <- calcNormFactors(y)

# Design matrix
design <- model.matrix(~0+group)
colnames(design) <- levels(group)

# Check design matrix rank and group sizes
design_rank <- qr(design)$rank
min_samples_per_group <- min(table(group))

# Use simplified mode if rank < columns or samples/group < 2
use_simplified_mode <- FALSE

if (design_rank < ncol(design) || min_samples_per_group < 2) {
  cat("Warning: fewer than 2 samples per group or insufficient rank\n")
  cat("Using simplified differential analysis mode, set BCV=0.2\n")
  use_simplified_mode <- TRUE

  # 切换到简化的exactTest模式而不是GLM模式
  bcv <- 0.2  # Set biological coefficient of variation
  y$common.dispersion <- bcv^2
  y$tagwise.dispersion <- rep(bcv^2, nrow(y))

  # 不调用glmQLFit，在差异分析时使用exactTest
} else {
  # 标准模式：正常估计离散度并进行似然比检验
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
}

# Identify control samples
control_samples <- sample_subset$sample_name[sample_subset$condition == "untreatment"]

# Identify treatment samples
treatment_groups <- sample_subset$sample_name[sample_subset$condition != "untreatment"]

# Compare each treatment against control
if (length(control_samples) > 0 && length(treatment_groups) > 0) {

  for (treatment in treatment_groups) {
    treatment_short <- gsub("-.*$", "", treatment)  # Shorten treatment name
    comparison_name <- paste0(treatment_short, "_vs_Control")

    # 获取该处理的分组标识
    treatment_group <- sample_subset$condition[sample_subset$sample_name == treatment]

    # 根据分析模式选择差异分析方法
    if (use_simplified_mode) {
      # 简化模式：使用exactTest
      # 注意：要使用exactTest，需要创建单组比较的设置
      # 编码处理组为1，对照组为2
      group_for_exact <- rep(2, ncol(y))  # Default all samples to control
      # Set 1 at treatment sample positions
      group_for_exact[y$samples$group == treatment_group] <- 1

      # 创建一个新的DGEList仅包含当前实验组和对照组
      samples_to_include <- y$samples$group %in% c(treatment_group, "untreatment")
      if (sum(samples_to_include) > 0) {
        y_exact <- y[, samples_to_include]
        y_exact$samples$group <- factor(ifelse(y_exact$samples$group == treatment_group, "treatment", "control"))

        # 使用exactTest进行差异分析
        et <- exactTest(y_exact, dispersion=bcv^2)
        results <- topTags(et, n = Inf)
      } else {
        stop("Error: unable to find ", treatment_group, " samples with untreatment")
      }
    } else {
      # 标准模式：使用glmQLFTest
      # 创建对比矩阵
      contrasts <- makeContrasts(contrasts = paste0(treatment_group, " - untreatment"), levels = design)

      # 差异分析
      qlf <- glmQLFTest(fit, contrast = contrasts)

      # 提取差异表达结果
      results <- topTags(qlf, n = Inf)
    }

    # 输出结果到文件
    result_file <- file.path(output_dir, paste0(comparison_name, "_results.tsv"))
    write.table(results$table, result_file, sep = "\t", quote = FALSE, row.names = TRUE)

    # 火山图
    results_df <- data.frame(results$table)
    results_df$TE <- rownames(results_df)
    results_df$significant <- ifelse(results_df$FDR < 0.05,
                                  ifelse(results_df$logFC > 0, "Up", "Down"), "Not Sig")

    p <- ggplot(results_df, aes(x = logFC, y = -log10(FDR), color = significant)) +
      geom_point(alpha = 0.7) +
      scale_color_manual(values = c("Down" = "blue", "Up" = "red", "Not Sig" = "gray")) +
      theme_bw() +
      labs(title = paste0(comparison_name, " differential analysis"),
           x = "log2 Fold Change",
           y = "-log10 FDR") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed")

    # 保存火山图
    ggsave(file.path(output_dir, paste0(comparison_name, "_volcano.pdf")), p, width = 8, height = 7)

    # 生成摘要
    sig_up <- sum(results_df$FDR < 0.05 & results_df$logFC > 0)
    sig_down <- sum(results_df$FDR < 0.05 & results_df$logFC < 0)

    cat(paste0("Differential expression analysis completed: ", comparison_name, "\n"))
    cat(paste0("Significantly upregulated TEs: ", sig_up, "\n"))
    cat(paste0("Significantly downregulated TEs: ", sig_down, "\n"))
  }
} else {
  cat("Error: no control or treatment samples found in sample config\n")
}

cat(paste0("All results saved to: ", output_dir, "\n"))

# Save session information
sessionInfo_file <- file.path(output_dir, "session_info.txt")
sink(sessionInfo_file)
sessionInfo()
sink()
