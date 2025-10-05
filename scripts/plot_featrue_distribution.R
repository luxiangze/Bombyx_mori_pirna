#!/usr/bin/env Rscript

# 检查并安装必要的包
packages <- c("tidyverse", "optparse", "tidyplots")
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)

# 加载必要的包
library(tidyverse)
library(optparse)
library(tidyplots)

# 解析命令行参数
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="输入CSV文件路径，第一列为类别，第二列为值", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="donut_plot.pdf",
              help="输出图片文件路径，支持PDF、PNG等格式 [默认= %default]", metavar="character")
)

# 创建选项解析器
opt_parser <- OptionParser(option_list=option_list,
                          description="\n使用CSV数据绘制环状图（Donut Plot）。\n输入文件格式: 第一列为类别名称，第二列为对应的数值。")

# 解析参数
opt <- parse_args(opt_parser)

# 检查输入文件是否提供
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("必须提供输入文件路径。使用 -i 或 --input 参数。", call.=FALSE)
}

# 检查输入文件是否存在
if (!file.exists(opt$input)) {
  stop(paste("输入文件不存在:", opt$input), call.=FALSE)
}

# 读取数据
data <- read_csv(opt$input, col_names = TRUE)

# 转换数据，使其符合tidyplot的要求
data_tidy <- data %>%
  pivot_longer(
    cols = -key1, # 除了 key1 的其他所有列
    names_to = "Group", # 将列名变成一个新变量
    values_to = "Expression" # 对应的值变成一个变量
  )

# 使用tidyplots创建环状图
data_tidy |>
  # 创建基础tidyplot对象，指定颜色方案
  tidyplot(y = Expression, color = key1) |>
  # 添加环状图，指定宽度和是否反转排序
  add_donut() |>
  adjust_size(width = 250, height = 250) |>
  split_plot(by = Group) |>
  # 保存图片
  save_plot(opt$output)

# 输出信息
cat(paste("\n环状图已保存到:", normalizePath(opt$output), "\n"))
cat(paste("包含类别数量:", nrow(data), "\n"))
