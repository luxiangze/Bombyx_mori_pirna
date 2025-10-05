#!/bin/bash

# 输入目录（包含fq.gz文件）
target_dir=$1

# 输出文件
OUTPUT_FILE="$target_dir/length_distribution.txt"

# 创建临时目录
TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

echo "Processing files in $target_dir..."

# 获取所有文件的列表
FILES=("$target_dir"/*fq.gz)
if [ ${#FILES[@]} -eq 0 ]; then
    echo "No fq.gz files found in $target_dir"
    exit 1
fi

# 清空或创建所有长度文件
> "$TMP_DIR/all_lengths.txt"

# 为每个文件创建一个临时文件来存储长度计数
for file in "${FILES[@]}"; do
    base_filename=$(basename "$file" .fq.gz)
    echo "Processing $base_filename..."
    
    # 创建一个空的计数文件
    > "$TMP_DIR/$base_filename.counts"
    
    # 提取每个序列的长度并计数
    zcat "$file" | awk 'NR % 4 == 2 {lengths[length($0)]++} 
                         END {for (l in lengths) print l, lengths[l]}' > "$TMP_DIR/$base_filename.tmp"
    
    # 将长度添加到总长度列表
    cut -d' ' -f1 "$TMP_DIR/$base_filename.tmp" >> "$TMP_DIR/all_lengths.txt"
    
    # 整理计数结果
    while read -r len count; do
        echo "$len $count" >> "$TMP_DIR/$base_filename.counts"
    done < "$TMP_DIR/$base_filename.tmp"
    
    rm "$TMP_DIR/$base_filename.tmp"
done

# 找出所有唯一长度并排序
sort -n -u "$TMP_DIR/all_lengths.txt" > "$TMP_DIR/unique_lengths.txt"

# 创建输出文件的表头
echo -n "Length" > "$OUTPUT_FILE"
for file in "${FILES[@]}"; do
    # 获取不带后缀的文件名
    base_filename=$(basename "$file" .fq.gz)
    # 确保移除 .fq.gz 后缀
    clean_name=${base_filename%_fq.gz}
    echo -n -e "\t$clean_name" >> "$OUTPUT_FILE"
done
echo "" >> "$OUTPUT_FILE"

# 填充表格
echo "Creating the distribution table..."
while read -r length; do
    echo -n "$length" >> "$OUTPUT_FILE"
    
    for file in "${FILES[@]}"; do
        base_filename=$(basename "$file" .fq.gz)
        # 查找此长度的计数
        count=$(grep -w "^$length" "$TMP_DIR/$base_filename.counts" | awk '{print $2}' || echo 0)
        # 如果没有找到计数，设为0
        if [ -z "$count" ]; then
            count=0
        fi
        echo -n -e "\t$count" >> "$OUTPUT_FILE"
    done
    
    echo "" >> "$OUTPUT_FILE"
done < "$TMP_DIR/unique_lengths.txt"

echo "Done! Length distribution saved to $OUTPUT_FILE"