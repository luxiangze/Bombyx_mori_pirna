#!/usr/bin/env Rscript

# 绘制piRNA表达散点图脚本
# 功能：
# 1. 计算各个piRNA的TPM值
# 2. 绘制散点图，横坐标为对照组，纵坐标为实验组
# 3. 根据p值和logFC设置点的颜色

# 检查并安装所需的R包
required_packages <- c("ggplot2", "dplyr", "optparse", "viridis", "scales")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# 加载所需的R包
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(optparse)
  library(viridis)
  library(scales)
})

# 定义命令行参数
option_list <- list(
  make_option(c("-d", "--data_dir"), type = "character", default = NULL,
              help = "包含edgeR结果文件的目录路径"),
  make_option(c("-o", "--output_dir"), type = "character", default = NULL,
              help = "输出图像的目录路径"),
  make_option(c("-p", "--p_value_cutoff"), type = "numeric", default = 0.01,
              help = "P值或FDR显著性阈值，默认为0.01"),
  make_option(c("--use_fdr"), action = "store_true", default = FALSE,
              help = "使用FDR而非P值作为显著性判断标准"),
  make_option(c("--min_length"), type = "numeric", default = 0,
              help = "序列最小长度过滤，默认为0（不过滤）"),
  make_option(c("--max_length"), type = "numeric", default = 1000,
              help = "序列最大长度过滤，默认为1000"),
  make_option(c("--pdf"), action = "store_true", default = FALSE,
              help = "输出PDF格式的图像，默认为PNG格式")
)

# 解析命令行参数
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# 检查必需的参数
if (is.null(opt$data_dir)) {
  stop("请提供数据目录路径 (-d 或 --data_dir)")
}

# 如果未指定输出目录，则使用当前目录
if (is.null(opt$output_dir)) {
  opt$output_dir <- getwd()
  message(paste("未指定输出目录，将使用当前目录:", opt$output_dir))
}

# 确保输出目录存在
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# 计算TPM值的函数
calculate_tpm <- function(counts) {
  # 确保counts是数值向量
  counts <- as.numeric(counts)
  
  # 计算每百万reads的比例因子
  scaling_factor <- sum(counts) / 1e6
  
  # 避免除以零
  if (scaling_factor == 0) {
    scaling_factor <- 1
  }
  
  # 计算TPM
  tpm <- counts / scaling_factor
  
  return(tpm)
}

# 处理文件并绘制散点图
process_file <- function(file_path, output_dir, p_value_cutoff, output_pdf = FALSE, use_fdr = FALSE, min_length = 0, max_length = 1000) {
  # 读取数据
  message(paste("处理文件:", basename(file_path)))
  data <- read.csv(file_path, row.names = 1, check.names = FALSE)
  
  # 计算序列长度并过滤
  seq_lengths <- nchar(rownames(data))
  message(paste("过滤前序列数量:", nrow(data)))
  data_filtered <- data[seq_lengths >= min_length & seq_lengths <= max_length, ]
  message(paste("过滤后序列数量:", nrow(data_filtered)))
  
  # 如果过滤后没有序列，返回错误
  if (nrow(data_filtered) == 0) {
    stop("过滤后没有序列符合长度要求")
  }
  
  # 使用过滤后的数据
  data <- data_filtered
  
  # 提取样本名称
  file_name <- basename(file_path)
  sample_names <- strsplit(gsub("_edger_results_with_counts.csv", "", file_name), "_vs_")[[1]]
  control_name <- sample_names[1]
  treatment_name <- sample_names[2]
  
  # 计算TPM值
  control_tpm <- calculate_tpm(data$control)
  treatment_tpm <- calculate_tpm(data$treatment)
  
  # 创建包含TPM值的数据框
  plot_data <- data.frame(
    piRNA = rownames(data),
    control_tpm = control_tpm,
    treatment_tpm = treatment_tpm,
    logFC = data$logFC,
    PValue = data$PValue,
    FDR = data$FDR
  )
  
  # 根据用户选择使用P值或FDR作为显著性判断标准
  if (use_fdr) {
    significance_value <- plot_data$FDR
    significance_type <- "FDR"
  } else {
    significance_value <- plot_data$PValue
    significance_type <- "P-value"
  }
  
  # 添加颜色分类
  plot_data$color_category <- ifelse(significance_value <= p_value_cutoff, 
                                     ifelse(plot_data$logFC > 0, "Up", 
                                            ifelse(plot_data$logFC < 0, "Down", "Neutral")),
                                     "NonSig")
  
  # 为了更好的可视化，对TPM值进行log10转换（添加小值避免log(0)）
  plot_data$log10_control_tpm <- log10(plot_data$control_tpm + 1)
  plot_data$log10_treatment_tpm <- log10(plot_data$treatment_tpm + 1)
  
  # 限制logFC范围为[-6, 6]
  plot_data$logFC_limited <- pmax(pmin(plot_data$logFC, 6), -6)
  
  # 计算上调和下调的点数
  up_count <- sum(plot_data$color_category == "Up")
  down_count <- sum(plot_data$color_category == "Down")
  
  # 创建散点图
  p <- ggplot(plot_data, aes(x = log10_control_tpm, y = log10_treatment_tpm)) +
    # 先绘制非显著性点（灰色）
    geom_point(data = subset(plot_data, color_category == "NonSig"), 
               color = "grey80", alpha = 0.4, size = 1) +
    # 再绘制显著性点（按logFC着色）
    geom_point(data = subset(plot_data, color_category != "NonSig"), 
               aes(color = logFC_limited), alpha = 0.4, size = 1.2) +
    # 添加对角线
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
    # 设置颜色映射，简化为蓝白红
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
    # 设置坐标轴标签
    labs(
      x = paste("log10(TPM) -", control_name),
      y = paste("log10(TPM) -", treatment_name),
      title = paste(treatment_name, "vs", control_name, "piRNA Expression"),
      subtitle = paste(significance_type, "cutoff:", p_value_cutoff)
    ) +
    # 添加上调和下调点的数量标注
    annotate("text", x = min(plot_data$log10_control_tpm, na.rm = TRUE) + 1, 
             y = max(plot_data$log10_treatment_tpm, na.rm = TRUE) - 1, 
             label = paste("Up: n =", up_count), 
             color = "red", size = 4, hjust = 0) +
    annotate("text", x = max(plot_data$log10_control_tpm, na.rm = TRUE) - 1, 
             y = min(plot_data$log10_treatment_tpm, na.rm = TRUE) + 1, 
             label = paste("Down: n =", down_count), 
             color = "blue", size = 4, hjust = 1) +
    # 设置主题
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  # 保存图像
  # 构建文件名，包含长度过滤信息
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
  
  message(paste("图像已保存至:", output_file))
  
  # 返回处理的数据
  return(list(plot = p, data = plot_data))
}

# 主函数
main <- function() {
  # 获取数据目录中所有的edgeR结果文件
  data_files <- list.files(opt$data_dir, 
                          pattern = "*_edger_results_with_counts.csv$", 
                          full.names = TRUE, 
                          recursive = FALSE)
  
  if (length(data_files) == 0) {
    stop(paste("在目录中未找到edgeR结果文件:", opt$data_dir))
  }
  
  message(paste("找到", length(data_files), "个文件进行处理"))
  
  # 处理每个文件
  results <- lapply(data_files, function(file) {
    tryCatch({
      process_file(file, opt$output_dir, opt$p_value_cutoff, opt$pdf, opt$use_fdr, opt$min_length, opt$max_length)
    }, error = function(e) {
      message(paste("处理文件时出错:", basename(file), "- 错误:", e$message))
      NULL
    })
  })
  
  # 过滤掉NULL结果
  results <- results[!sapply(results, is.null)]
  
  if (length(results) > 0) {
    message(paste("成功处理了", length(results), "个文件"))
  } else {
    message("没有成功处理任何文件")
  }
}

# 执行主函数
main()
