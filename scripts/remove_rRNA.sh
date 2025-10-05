#!/bin/bash

# 脚本：remove_rRNA.sh
# 功能：从小RNA数据中移除rRNA序列
# 用法：bash remove_rRNA.sh <输入目录> <bowtie索引目录> <输出目录>

# 检查参数数量
if [ $# -ne 3 ]; then
    echo "错误：参数数量不正确"
    echo "用法：bash $0 <输入目录> <bowtie索引目录> <输出目录>"
    echo "示例：bash $0 /path/to/input /path/to/bowtie_index /path/to/output"
    exit 1
fi

# 获取参数
INPUT_DIR="$1"
BOWTIE_INDEX_DIR="$2"
OUTPUT_DIR="$3"

# 创建日志目录
PROJECT_DIR=$(dirname $(dirname "$0"))
LOG_DIR="$PROJECT_DIR/logs"
mkdir -p "$LOG_DIR"

# 创建日志文件（使用时间戳）
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="$LOG_DIR/remove_rRNA_${TIMESTAMP}.log"

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

# 记录脚本开始执行的信息
echo "开始执行rRNA移除任务：$(date)" | tee -a "$LOG_FILE"
echo "输入目录：$INPUT_DIR" | tee -a "$LOG_FILE"
echo "Bowtie索引目录：$BOWTIE_INDEX_DIR" | tee -a "$LOG_FILE"
echo "输出目录：$OUTPUT_DIR" | tee -a "$LOG_FILE"
echo "-----------------------------------" | tee -a "$LOG_FILE"

# 检查输入目录是否存在
if [ ! -d "$INPUT_DIR" ]; then
    echo "错误：输入目录 '$INPUT_DIR' 不存在" | tee -a "$LOG_FILE"
    exit 1
fi

# 检查bowtie索引目录是否存在
if [ ! -d "$BOWTIE_INDEX_DIR" ]; then
    echo "错误：Bowtie索引目录 '$BOWTIE_INDEX_DIR' 不存在" | tee -a "$LOG_FILE"
    exit 1
fi

# 查找所有输入文件（包括gz压缩文件）
echo "查找输入文件..." | tee -a "$LOG_FILE"
FILE_LIST=$(mktemp)
find "$INPUT_DIR" -type f \( -name "*.fa" -o -name "*.fasta" -o -name "*.fastq" -o -name "*.fq" -o -name "*.fa.gz" -o -name "*.fasta.gz" -o -name "*.fastq.gz" -o -name "*.fq.gz" \) > "$FILE_LIST"

# 统计文件数量
FILE_COUNT=$(wc -l < "$FILE_LIST")
echo "找到 $FILE_COUNT 个文件需要处理" | tee -a "$LOG_FILE"

# 获取CPU核心数（使用90%的可用核心）
TOTAL_CORES=$(nproc)
USABLE_CORES=$((TOTAL_CORES * 9 / 10))
if [ $USABLE_CORES -lt 1 ]; then
    USABLE_CORES=1
fi

# 根据文件数量和核心数量计算每个文件分配的核心数
if [ $FILE_COUNT -gt 0 ]; then
    # 最大并行任务数，默认为文件数量，但不超过可用核心数
    MAX_PARALLEL_JOBS=$FILE_COUNT
    if [ $MAX_PARALLEL_JOBS -gt $USABLE_CORES ]; then
        MAX_PARALLEL_JOBS=$USABLE_CORES
    fi
    
    # 每个文件分配的核心数
    CORES_PER_FILE=1
    if [ $FILE_COUNT -lt $USABLE_CORES ]; then
        CORES_PER_FILE=$((USABLE_CORES / FILE_COUNT))
    fi
    
    echo "并行处理文件数：$MAX_PARALLEL_JOBS，每个文件使用核心数：$CORES_PER_FILE" | tee -a "$LOG_FILE"
else
    echo "警告：未找到文件可处理" | tee -a "$LOG_FILE"
    rm -f "$FILE_LIST"
    exit 0
fi

# 定义处理单个文件的函数
process_file() {
    local file=$1
    local cores=$2
    local log_file=$3
    
    # 获取文件名和扩展名
    local filename=$(basename "$file")
    local is_gzipped=false
    local extension=""
    local basename=""
    
    # 检查是否为gz压缩文件
    if [[ "$filename" == *.gz ]]; then
        is_gzipped=true
        # 移除.gz扩展名获取原始文件名
        local orig_filename=${filename%.gz}
        extension="${orig_filename##*.}"
        basename="${orig_filename%.*}"
    else
        extension="${filename##*.}"
        basename="${filename%.*}"
    fi
    
    # 输出文件路径
    local output_file="$OUTPUT_DIR/${basename}.${extension}"
    
    # 如果原文件是压缩的，则输出也压缩
    if [ "$is_gzipped" = true ]; then
        output_file="${output_file}.gz"
    fi
    
    # 创建临时文件用于输出
    local temp_output="$(mktemp -p /tmp "${basename}_x_rRNA.XXXXXX")"
    
    # 记录处理信息
    echo "[进程 $BASHPID] 开始处理文件：$filename" >> "$log_file"
    
    # 根据文件扩展名确定输入格式
    local format_option=""
    if [[ "$extension" == "fa" || "$extension" == "fasta" ]]; then
        format_option="-f"
    elif [[ "$extension" == "fq" || "$extension" == "fastq" ]]; then
        format_option="-q"
    else
        echo "[进程 $BASHPID]   警告：无法确定文件格式，默认使用FASTA格式" >> "$log_file"
        format_option="-f"
    fi
    
    # Bowtie支持直接读取gz压缩文件，不需要解压
    local input_for_bowtie="$file"
    
    # 运行bowtie，提取未比对到rRNA的序列
    echo "[进程 $BASHPID]   运行Bowtie比对..." >> "$log_file"
    bowtie -p $cores --un "$temp_output" $format_option \
           -v 2 -a --best --strata --norc \
           -S -x "$BOWTIE_INDEX_DIR/rRNA" \
           "$input_for_bowtie" \
           /dev/null 2>> "$log_file"
    
    # 检查bowtie运行状态
    if [ $? -eq 0 ]; then
        # 如果需要压缩输出
        if [ "$is_gzipped" = true ]; then
            echo "[进程 $BASHPID]   压缩输出文件..." >> "$log_file"
            gzip -c "$temp_output" > "$output_file"
        else
            mv "$temp_output" "$output_file"
        fi
        echo "[进程 $BASHPID]   成功：已将非rRNA序列保存到 $output_file" >> "$log_file"
    else
        echo "[进程 $BASHPID]   错误：处理文件 $file 时出错" >> "$log_file"
    fi
    
    # 清理临时文件
    if [ -f "$temp_output" ]; then
        rm -f "$temp_output"
    fi
    
    echo "[进程 $BASHPID] 完成处理文件：$filename" >> "$log_file"
    echo "[进程 $BASHPID] -----------------------------------" >> "$log_file"
}

# 开始并行处理文件
echo "开始处理文件..." | tee -a "$LOG_FILE"

# 初始化计数器
COMPLETED=0

# 读取文件列表并处理
while read -r file; do
    # 检查当前运行的后台任务数
    while [ $(jobs -p | wc -l) -ge $MAX_PARALLEL_JOBS ]; do
        # 等待任一后台任务完成
        sleep 1
    done
    
    # 获取文件名
    filename=$(basename "$file")
    echo "启动处理文件：$filename (进度: $COMPLETED/$FILE_COUNT)" | tee -a "$LOG_FILE"
    
    # 在后台运行处理函数
    process_file "$file" $CORES_PER_FILE "$LOG_FILE" &
    
    # 增加计数器
    COMPLETED=$((COMPLETED + 1))
done < "$FILE_LIST"

# 等待所有后台任务完成
echo "等待所有任务完成..." | tee -a "$LOG_FILE"
wait

# 清理临时文件列表
rm -f "$FILE_LIST"

echo "rRNA移除任务完成：$(date)" | tee -a "$LOG_FILE"
echo "结果已保存到：$OUTPUT_DIR" | tee -a "$LOG_FILE"
echo "日志文件：$LOG_FILE" | tee -a "$LOG_FILE"

exit 0
