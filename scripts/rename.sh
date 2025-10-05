#!/bin/bash

# 使用方法检查
if [ $# -ne 2 ]; then
    echo "用法: $0 <文件目录> <样本映射文件>"
    echo "例如: $0 data/raw_data data/sample_map.txt"
    exit 1
fi

# 获取参数
DATA_DIR=$1
SAMPLE_MAP=$2

# 检查目录和映射文件是否存在
if [ ! -d "$DATA_DIR" ]; then
    echo "错误: 目录 '$DATA_DIR' 不存在"
    exit 1
fi

if [ ! -f "$SAMPLE_MAP" ]; then
    echo "错误: 样本映射文件 '$SAMPLE_MAP' 不存在"
    exit 1
fi

# 创建md5文件
MD5_FILE="$DATA_DIR/md5.txt"
> "$MD5_FILE"

# 样本名称映射文件处理
while read -r old_name new_name; do
  # 查找原始文件名的匹配文件
  for file in "$DATA_DIR/${old_name}"*.fastq.gz; do
    if [[ -f "$file" ]]; then
      # 重命名文件
      new_file="$DATA_DIR/${new_name}.fq.gz"
      echo "重命名: $file -> $new_file"
      mv "$file" "$new_file"
      
      # 计算新的 md5 值
      md5sum "$new_file" >> "$MD5_FILE"
    fi
  done
done < "$SAMPLE_MAP"

echo "所有文件已重命名，md5 校验值已更新到 $MD5_FILE"