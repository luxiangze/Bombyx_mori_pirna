#!/usr/bin/env Rscript

# 用法: Rscript pirna_edger_pairwise.R control_file.uniq.fa treatment_file.uniq.fa

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 2){
  cat("用法: Rscript pirna_edger_pairwise.R control_file.uniq.fa treatment_file.uniq.fa\n")
  quit(status=1)
}

control_file <- args[1]
treatment_file <- args[2]

# 加载必要的包
if (!require(edgeR, quietly=TRUE)) {
  install.packages("edgeR")
  library(edgeR)
}

# 安装并加载Biostrings包（用于读取FASTA文件）
if (!require(Biostrings, quietly=TRUE)) {
  if (!require(BiocManager, quietly=TRUE))
    install.packages("BiocManager")
  BiocManager::install("Biostrings", update=FALSE, ask=FALSE)
  library(Biostrings)
}

# 读取uniq.fa，返回一个数据框（geneid=序列，count=reads数）
read_fa_counts <- function(fafile) {
  # 使用Biostrings读取FASTA文件
  fa <- readDNAStringSet(fafile)
  
  # 提取序列和ID
  sequences <- as.character(fa)
  headers <- names(fa)
  
  # 从header中提取reads数
  counts <- integer(length(headers))
  for (i in seq_along(headers)) {
    parts <- unlist(strsplit(headers[i], "-"))
    counts[i] <- as.integer(parts[2])
  }
  
  # 创建命名向量，以序列为名称
  names(counts) <- sequences
  return(counts)
}

# 读取两个样本
control_counts <- read_fa_counts(control_file)
treat_counts <- read_fa_counts(treatment_file)

# 合并所有geneid
all_geneids <- union(names(control_counts), names(treat_counts))
control_vec <- setNames(rep(0, length(all_geneids)), all_geneids)
treat_vec <- control_vec
control_vec[names(control_counts)] <- control_counts
treat_vec[names(treat_counts)] <- treat_counts

count_matrix <- cbind(control=control_vec, treatment=treat_vec)
rownames(count_matrix) <- all_geneids

# 分组信息
group <- factor(c("control", "treatment"), levels=c("control", "treatment"))

# edgeR分析
dge <- DGEList(counts = count_matrix, group = group)
dge <- calcNormFactors(dge)

# 由于没有重复，手动设置离散度（dispersion）
# 使用常见的无重复数据分析方法
cat("注意：数据没有生物学重复，使用固定离散度0.1进行分析\n")
cat("结果仅供参考，建议进行实验验证\n")

# 设置一个合理的离散度值
dge$common.dispersion <- 0.1

# 执行精确检验
et <- exactTest(dge, pair=c("control", "treatment"))
res <- topTags(et, n=Inf)$table

# 输出结果
# 创建输出文件名（基于输入文件名）
control_name <- sub("\\.uniq\\.fa$", "", basename(control_file))
treatment_name <- sub("\\.uniq\\.fa$", "", basename(treatment_file))
output_prefix <- paste0("results/pirna_edger_pairwise_out/", control_name, "_vs_", treatment_name)

# 保存表达矩阵和差异分析结果
count_file <- paste0(output_prefix, "_count_matrix.csv")
result_file <- paste0(output_prefix, "_edger_results.csv")

write.csv(count_matrix, file=count_file, quote=F)
write.csv(res, file=result_file, quote=F)

# 输出简单统计信息
cat("分析完成！\n")
cat("总序列数：", nrow(count_matrix), "\n")
cat("显著差异序列数 (FDR < 0.05)：", sum(res$FDR < 0.05), "\n")
cat("上调序列数 (logFC > 0)：", sum(res$FDR < 0.05 & res$logFC > 0), "\n")
cat("下调序列数 (logFC < 0)：", sum(res$FDR < 0.05 & res$logFC < 0), "\n")
cat("结果已保存为：\n")
cat("  - 表达矩阵：", count_file, "\n")
cat("  - 差异分析结果：", result_file, "\n")