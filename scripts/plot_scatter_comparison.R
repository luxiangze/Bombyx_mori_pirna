#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
#
# 通用散点图比较工具
# 用于绘制两个样本之间的散点图比较，支持多种自定义选项
#
# 使用方法:
# Rscript plot_scatter_comparison.R [参数]
#
# 参数:
#   --input_dir=<目录>       输入CSV文件所在目录
#   --output_dir=<目录>      输出PDF文件保存目录
#   --pattern=<模式>         用于筛选CSV文件的模式，默认为"_pirna_trim.csv"
#   --id_column=<列名>       用于合并数据的ID列名，默认为"piRNA_ID"
#   --value_column=<列名>    用于比较的数值列名，默认为"Trim_index"
#   --log_scale=<TRUE/FALSE> 是否使用对数坐标轴，默认为TRUE
#   --point_size=<数值>      点的大小，默认为1
#   --point_alpha=<数值>     点的透明度，默认为0.7
#   --width=<数值>           输出PDF的宽度(英寸)，默认为7
#   --height=<数值>          输出PDF的高度(英寸)，默认为7
#   --help                   显示帮助信息

# 检查并安装缺失的包
packages <- c("ggplot2", "dplyr", "tidyr", "readr", "optparse", "scales", "viridis")
missing_packages <- packages[!packages %in% installed.packages()[,"Package"]]
if(length(missing_packages) > 0) {
  message(paste0("正在安装缺失的包: ", paste(missing_packages, collapse=", ")))
  install.packages(missing_packages, repos="https://cloud.r-project.org")
}

# 加载必要的库
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(optparse)
  library(scales)
  library(viridis)
})

# 定义命令行参数
option_list <- list(
  make_option("--input_dir", type="character", default=NULL,
              help="输入CSV文件所在目录"),
  make_option("--output_dir", type="character", default=NULL,
              help="输出PDF文件保存目录"),
  make_option("--pattern", type="character", default="_pirna_trim.csv",
              help="用于筛选CSV文件的模式 [默认: %default]"),
  make_option("--id_column", type="character", default="piRNA_ID",
              help="用于合并数据的ID列名 [默认: %default]"),
  make_option("--value_column", type="character", default="Trim_index",
              help="用于比较的数值列名 [默认: %default]"),
  make_option("--log_scale", type="logical", default=TRUE,
              help="是否使用对数坐标轴 [默认: %default]"),
  make_option("--point_size", type="double", default=1,
              help="点的大小 [默认: %default]"),
  make_option("--point_alpha", type="double", default=0.7,
              help="点的透明度 [默认: %default]"),
  make_option("--width", type="double", default=7,
              help="输出PDF的宽度(英寸) [默认: %default]"),
  make_option("--height", type="double", default=7,
              help="输出PDF的高度(英寸) [默认: %default]"),
  make_option("--sample_config", type="character", default=NULL, 
              help="样本配置文件路径，包含样本名称和类型信息")
)

# 解析命令行参数
opt_parser <- OptionParser(option_list=option_list,
                          description="通用散点图比较工具 - 用于绘制两个样本之间的散点图比较")
opt <- parse_args(opt_parser)

# 检查必要的参数
if (is.null(opt$input_dir) || is.null(opt$output_dir)) {
  print_help(opt_parser)
  stop("必须指定输入和输出目录", call.=FALSE)
}

# 创建输出目录
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# 简化样本名称的函数
simplify_sample_name <- function(sample_name) {
  # 去除常见的后缀
  suffixes <- c(".map", ".fa.rmdup.map", ".fa.collapsed.map", ".fa.collapsed.no-dust.map",
                ".collapsed.no-dust.map", ".no-dust.map", ".fa.collapsed.no-dust",
                ".collapsed.no-dust", ".no-dust", ".fa")

  for (suffix in suffixes) {
    if (endsWith(sample_name, suffix)) {
      sample_name <- substr(sample_name, 1, nchar(sample_name) - nchar(suffix))
      break
    }
  }

  # 移除测序相关的后缀
  patterns <- c("_1.fq", "_2.fq", ".fq", "_1.fastq", "_2.fastq", ".fastq")
  for (pattern in patterns) {
    sample_name <- gsub(pattern, "", sample_name)
  }

  # 如果样本名称包含"-KD"，只保留该部分
  if (grepl("-KD", sample_name)) {
    parts <- strsplit(sample_name, "-KD")[[1]]
    if (length(parts) > 1) {
      kd_parts <- strsplit(parts[2], "_")[[1]]
      sample_name <- paste0(parts[1], "-KD", kd_parts[1])
    }
  }

  return(sample_name)
}

# 获取输入文件列表
input_files <- list.files(path = opt$input_dir, pattern = opt$pattern, full.names = TRUE)
if (length(input_files) < 2) {
  stop("至少需要2个样本才能进行比较分析", call.=FALSE)
}

# 读取样本配置文件
sample_info <- list()
if (!is.null(opt$sample_config) && file.exists(opt$sample_config)) {
  message(paste("读取样本配置文件:", opt$sample_config))
  config_lines <- readLines(opt$sample_config)
  for (line in config_lines) {
    if (line != "") {
      parts <- strsplit(line, " ")[[1]]
      if (length(parts) >= 2) {
        sample_name <- parts[1]
        sample_type <- parts[2]
        sample_info[[sample_name]] <- sample_type
        message(paste("  样本:", sample_name, ", 类型:", sample_type))
      }
    }
  }
  message(paste("样本配置文件中包含", length(sample_info), "个样本"))
}

# 读取所有样本数据
all_data <- list()
for (file in input_files) {
  # 从文件名中提取样本名
  file_name <- basename(file)
  sample_name <- gsub(opt$pattern, "", file_name)

  # 简化样本名称
  simple_name <- simplify_sample_name(sample_name)

  # 读取CSV文件
  message(paste("读取文件:", file_name))
  df <- read_csv(file, show_col_types = FALSE)

  # 检查必要的列是否存在
  if (!opt$id_column %in% colnames(df) || !opt$value_column %in% colnames(df)) {
    warning(paste("文件", file_name, "缺少必要的列，跳过"))
    next
  }

  # 存储数据
  all_data[[simple_name]] <- df
  message(paste("  行数:", nrow(df)))
  message(paste("  列名:", paste(colnames(df), collapse=", ")))
}

message(paste("\n总共找到", length(all_data), "个样本"))

# 绘制每对样本的比较图
message("\
开始绘制样本比较图...")
sample_names <- names(all_data)
n_samples <- length(sample_names)

# 定义一个函数来获取样本类型
get_sample_type <- function(sample_name) {
  # 首先尝试直接匹配
  if (sample_name %in% names(sample_info)) {
    return(sample_info[[sample_name]])
  }

  # 如果没有直接匹配，尝试部分匹配
  for (name in names(sample_info)) {
    if (grepl(name, sample_name, fixed = TRUE) || grepl(sample_name, name, fixed = TRUE)) {
      return(sample_info[[name]])
    }
  }

  # 如果没有匹配，返回默认值
  return("unknown")
}

# 定义一个函数来决定样本的顺序
# 对照组应该在x轴，处理组应该在y轴
should_swap_samples <- function(sample1, sample2) {
  type1 <- get_sample_type(sample1)
  type2 <- get_sample_type(sample2)

  # 如果第一个样本是处理组，第二个是对照组，则交换
  if (type1 == "treatment" && type2 == "untreatment") {
    return(TRUE)
  }

  # 否则不交换
  return(FALSE)
}

# 定义一个函数来检查是否是对照组和处理组的比较
is_control_vs_treatment <- function(sample1, sample2) {
  type1 <- get_sample_type(sample1)
  type2 <- get_sample_type(sample2)

  # 如果一个是对照组，一个是处理组，则返回TRUE
  if ((type1 == "untreatment" && type2 == "treatment") ||
      (type1 == "treatment" && type2 == "untreatment")) {
    return(TRUE)
  }

  # 否则返回FALSE
  return(FALSE)
}

for (i in 1:(n_samples-1)) {
  for (j in (i+1):n_samples) {
    sample1_name <- sample_names[i]
    sample2_name <- sample_names[j]

    # 检查是否是对照组和处理组的比较
    if (!is_control_vs_treatment(sample1_name, sample2_name)) {
      message(paste("  跳过:", sample1_name, "vs", sample2_name, "(非对照组与处理组的比较)"))
      next
    }

    # 检查是否需要交换样本顺序
    if (should_swap_samples(sample1_name, sample2_name)) {
      temp <- sample1_name
      sample1_name <- sample2_name
      sample2_name <- temp
    }

    message(paste("  绘制:", sample1_name, "vs", sample2_name))

    # 获取样本类型
    type1 <- get_sample_type(sample1_name)
    type2 <- get_sample_type(sample2_name)

    # 获取两个样本的数据
    df1 <- all_data[[sample1_name]]
    df2 <- all_data[[sample2_name]]

    # 合并两个样本的数据
    merged_df <- df1 %>%
      select(!!sym(opt$id_column), !!sym(opt$value_column)) %>%
      inner_join(
        df2 %>% select(!!sym(opt$id_column), !!sym(opt$value_column)),
        by = opt$id_column,
        suffix = c("_1", "_2")
      )

    # 计算fold change
    merged_df <- merged_df %>%
      mutate(fold_change = .data[[paste0(opt$value_column, "_2")]] / (.data[[paste0(opt$value_column, "_1")]] + 1e-10))

    # 计算fold change的范围，用于颜色映射
    fc_min <- min(merged_df$fold_change)
    fc_max <- max(merged_df$fold_change)

    # 确保范围是对称的（相对于1）
    if (fc_min < 1 && 1/fc_min > fc_max) {
      vmax <- 1/fc_min
      vmin <- fc_min
    } else {
      vmax <- fc_max
      vmin <- 1/fc_max
    }

    # 不再限制颜色映射范围
    # 记录fold change的范围供参考
    message(paste("    Fold change范围:", round(fc_min, 2), "-", round(fc_max, 2)))

    # 创建颜色映射函数 - 深蓝黄深红配色
    color_scale <- scale_color_gradient2(
      low = "darkblue",
      mid = "yellow",
      high = "darkred",
      midpoint = 1,
      trans = "log10",
      name = "Fold change"
    )

    # 获取数据范围
    min_val <- min(c(merged_df[[paste0(opt$value_column, "_1")]], merged_df[[paste0(opt$value_column, "_2")]]))
    max_val <- max(c(merged_df[[paste0(opt$value_column, "_1")]], merged_df[[paste0(opt$value_column, "_2")]]))

    # 确保最小值为正数，用于对数坐标轴
    min_val <- max(min_val, 1e-10)

    # 创建散点图
    p <- ggplot(merged_df, aes(
      x = .data[[paste0(opt$value_column, "_1")]],
      y = .data[[paste0(opt$value_column, "_2")]],
      color = fold_change
    )) +
      geom_point(size = opt$point_size, alpha = 0.4) + # 调整点的透明度为0.4
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

    # 根据需要设置对数坐标轴
    if (opt$log_scale) {
      p <- p +
        scale_x_log10(
          limits = c(min_val, max_val * 1.1),
          labels = function(x) format(x, scientific = FALSE, digits = 1) # 使用正确的对数刻度格式
        ) +
        scale_y_log10(
          limits = c(min_val, max_val * 1.1),
          labels = function(x) format(x, scientific = FALSE, digits = 1) # 使用正确的对数刻度格式
        )
    } else {
      p <- p +
        scale_x_continuous(limits = c(0, max_val * 1.1)) +
        scale_y_continuous(limits = c(0, max_val * 1.1))
    }

    # 不再添加样本类型的注释

    # 保存图表
    output_file <- file.path(opt$output_dir, paste0(sample1_name, "_vs_", sample2_name, ".pdf"))
    ggsave(output_file, plot = p, width = opt$width, height = opt$height)
    message(paste("    已保存到:", output_file))
  }
}

message(paste("\n绘图完成，结果保存在:", opt$output_dir))
