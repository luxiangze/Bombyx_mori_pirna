#!/usr/bin/env Rscript

# TE差异表达分析脚本 (使用edgeR)
# 用法: Rscript te_differential_analysis.R [计数矩阵文件] [样本配置文件] [输出目录]

# 加载必要的库
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

# 解析命令行参数
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

# 创建输出目录
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 读取样本配置信息
sample_info <- read.table(sample_config, header = FALSE, stringsAsFactors = FALSE)
colnames(sample_info) <- c("sample_name", "condition")

# 读取表达计数数据
counts_data <- read.delim(counts_file, header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

# 根据实际数据调整样本名映射
# 注意：这里根据实际数据文件中的列名进行硬编码映射
# 如果有多个重复样本，可能需要调整此映射
sample_mapping <- c(
  "control" = "Control",
  "DHX15" = "DHX15-KD3",
  "RNPS1" = "RNPS1-KD1",
  "SUGP1-Rep1" = "SUGP1-KD2",
  "SUGP1-Rep2" = "SUGP1-KD2"  # 注意Rep2被映射到同一个样本名
)

# 调整counts_data的列名，使其与sample_config中的名称匹配
colnames_original <- colnames(counts_data)
colnames_mapped <- sapply(colnames_original, function(x) {
  if (x %in% names(sample_mapping)) {
    return(sample_mapping[x])
  } else {
    return(x)
  }
})
colnames(counts_data) <- colnames_mapped

# 筛选配置文件中的样本
sample_subset <- subset(sample_info, sample_name %in% colnames_mapped)

# 如果需要，处理重复样本（例如对多个重复样本取平均值）
if (length(unique(colnames_mapped)) < length(colnames_mapped)) {
  cat("检测到重复样本，合并处理...\n")
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

# 确保counts_data的列顺序与sample_subset匹配
common_samples <- intersect(sample_subset$sample_name, colnames(counts_data))
counts_data <- counts_data[, common_samples]
sample_subset <- sample_subset[match(common_samples, sample_subset$sample_name), ]

# 创建DGEList对象
group <- factor(sample_subset$condition)
y <- DGEList(counts = counts_data, group = group)

# 过滤低表达的转座子
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]

# 标准化
y <- calcNormFactors(y)

# 设计矩阵
design <- model.matrix(~0+group)
colnames(design) <- levels(group)

# 判断规划矩阵是否有充分的自由度
design_rank <- qr(design)$rank
min_samples_per_group <- min(table(group))

# 如果设计矩阵的秩小于列数或每组样本数少于2，则使用简化模式
use_simplified_mode <- FALSE

if (design_rank < ncol(design) || min_samples_per_group < 2) {
  cat("警告：对照组或处理组样本数量少于2或计算秩不足\n")
  cat("使用简化的差异分析模式，设置BCV=0.2\n")
  use_simplified_mode <- TRUE

  # 切换到简化的exactTest模式而不是GLM模式
  bcv <- 0.2  # 设置生物学变异系数
  y$common.dispersion <- bcv^2
  y$tagwise.dispersion <- rep(bcv^2, nrow(y))

  # 不调用glmQLFit，在差异分析时使用exactTest
} else {
  # 标准模式：正常估计离散度并进行似然比检验
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
}

# 找出对照组样本
control_samples <- sample_subset$sample_name[sample_subset$condition == "untreatment"]

# 找出所有非对照组样本
treatment_groups <- sample_subset$sample_name[sample_subset$condition != "untreatment"]

# 对每个处理组与对照组进行比较
if (length(control_samples) > 0 && length(treatment_groups) > 0) {

  for (treatment in treatment_groups) {
    treatment_short <- gsub("-.*$", "", treatment)  # 简化处理组名称
    comparison_name <- paste0(treatment_short, "_vs_Control")

    # 获取该处理的分组标识
    treatment_group <- sample_subset$condition[sample_subset$sample_name == treatment]

    # 根据分析模式选择差异分析方法
    if (use_simplified_mode) {
      # 简化模式：使用exactTest
      # 注意：要使用exactTest，需要创建单组比较的设置
      # 编码处理组为1，对照组为2
      group_for_exact <- rep(2, ncol(y))  # 默认所有样本为对照组
      # 在处理组样本的位置设置为1
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
        stop("错误：无法找到", treatment_group, "与 untreatment 的样本")
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
      labs(title = paste0(comparison_name, " 差异分析"),
           x = "log2 Fold Change",
           y = "-log10 FDR") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed")

    # 保存火山图
    ggsave(file.path(output_dir, paste0(comparison_name, "_volcano.pdf")), p, width = 8, height = 7)

    # 生成摘要
    sig_up <- sum(results_df$FDR < 0.05 & results_df$logFC > 0)
    sig_down <- sum(results_df$FDR < 0.05 & results_df$logFC < 0)

    cat(paste0("差异表达分析完成: ", comparison_name, "\n"))
    cat(paste0("显著上调的转座子: ", sig_up, "\n"))
    cat(paste0("显著下调的转座子: ", sig_down, "\n"))
  }
} else {
  cat("错误：在样本配置中没有找到对照组或处理组样本\n")
}

cat(paste0("所有分析结果保存在: ", output_dir, "\n"))

# 保存会话信息
sessionInfo_file <- file.path(output_dir, "session_info.txt")
sink(sessionInfo_file)
sessionInfo()
sink()
