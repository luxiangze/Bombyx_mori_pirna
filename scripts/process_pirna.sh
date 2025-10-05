#!/bin/bash

# 自动处理piRNA数据的脚本
# 用法: scripts/process_pirna.sh <输入目录> <输出目录> <参考基因组路径> [--skip-existing] [--reverse]

################################################################################
# 函数设置
################################################################################

# 日志函数
log() {
    local message="[$(date '+%Y-%m-%d %H:%M:%S')] $1"
    echo "$message" | tee -a "$MAIN_LOG"
}

# 简化文件名函数
simplify_filename() {
    local filename=$1
    # 移除常见的后缀
    filename=$(echo "$filename" | sed -E 's/_1\.fq|_2\.fq|\.fq|_1\.fastq|_2\.fastq|\.fastq//g')
    echo "$filename"
}

# 处理单个文件的函数
process_file() {
    local FASTQ_FILE=$1
    local CORES_TO_USE=6  # 分配给该任务的CPU核心数
    local RAW_FILENAME=$(basename "$FASTQ_FILE" | sed 's/\.[^.]*$//')
    # 简化文件名
    local FILENAME=$(simplify_filename "$RAW_FILENAME")
    local FILE_LOG="$LOG_SUBDIR/${FILENAME}.log"
    
    # 检查最终输出文件是否已存在 (使用原始文件名检查，以避免冲突)
    local FINAL_OUTPUT="$OUTPUT_DIR/${FILENAME}.fa.collapsed.map"
    if $SKIP_EXISTING && [ -f "$FINAL_OUTPUT" ]; then
        log "跳过已处理的文件: $FILENAME (输出文件已存在)"
        return 0
    fi
    
    # 记录开始处理
    log "开始处理文件: $FILENAME (PID: $$)"
    
    # 1. 转换fastq到fasta
    log "[$FILENAME] 步骤1: 序列预处理..." >> "$FILE_LOG" 2>&1
    log "[$FILENAME] 使用 $CORES_TO_USE 个CPU核心运行seqkit" >> "$FILE_LOG" 2>&1
    
    # 根据是否需要反向处理选择不同的处理流程
    if $REVERSE_SEQ; then
        log "[$FILENAME] 启用序列反向处理" >> "$FILE_LOG" 2>&1
        # 先进行长度筛选，然后反向处理序列，最后进行序列折叠
        seqkit fq2fa "$FASTQ_FILE" -j $CORES_TO_USE | \
        seqkit seq -m 23 -j $CORES_TO_USE | \
        seqkit seq -r -j $CORES_TO_USE | \
        fastx_collapser -o "$OUTPUT_DIR/${FILENAME}.fa.collapsed" >> "$FILE_LOG" 2>&1
        
        # 使用已生成的反向参考基因组
        local REF_DIR=$(dirname "$PIRNA_REF")
        local REVERSE_REF="$REF_DIR/$(basename "$PIRNA_REF").reverse"
        
        # 2. piRNA映射（使用反向参考基因组）
        log "[$FILENAME] 步骤2: 进行piRNA映射（使用反向序列）..." >> "$FILE_LOG" 2>&1
        # 使用seqmap进行piRNA映射
        seqmap 0 "$OUTPUT_DIR/${FILENAME}.fa.collapsed" "$REVERSE_REF" "$OUTPUT_DIR/${FILENAME}.fa.collapsed.map" \
        /do_not_output_probe_without_match /forward_strand /cut:3,23 /output_all_matches >> "$FILE_LOG" 2>&1
        # 使用awk过滤出map文件的第1，2，4列
        awk '{print $4"\t"$2"\t"$1}' "$OUTPUT_DIR/${FILENAME}.fa.collapsed.map" > "$OUTPUT_DIR/${FILENAME}.fa.collapsed.map.filtered"
        mv "$OUTPUT_DIR/${FILENAME}.fa.collapsed.map.filtered" "$OUTPUT_DIR/${FILENAME}.fa.collapsed.map"
    else
        # 原始处理流程（不进行反向处理）
        seqkit fq2fa "$FASTQ_FILE" -j $CORES_TO_USE | \
        seqkit seq -m 26 -j $CORES_TO_USE | \
        fastx_collapser -o "$OUTPUT_DIR/${FILENAME}.fa.collapsed" >> "$FILE_LOG" 2>&1
        
        # 2. piRNA映射
        log "[$FILENAME] 步骤2: 进行piRNA映射..." >> "$FILE_LOG" 2>&1
        # 使用seqmap进行piRNA映射
        seqmap 0 "$PIRNA_REF" "$OUTPUT_DIR/${FILENAME}.fa.collapsed" "$OUTPUT_DIR/${FILENAME}.fa.collapsed.map" \
        /do_not_output_probe_without_match /forward_strand /cut:3,23 /output_all_matches >> "$FILE_LOG" 2>&1
        # 使用awk过滤出map文件的第1，2，4列，3减去第二列，然后过滤掉第二列小于-2的行
        awk '{print $4"\t"3-$2"\t"$1}' "$OUTPUT_DIR/${FILENAME}.fa.collapsed.map" | awk '$2 > -3' > "$OUTPUT_DIR/${FILENAME}.fa.collapsed.map.filtered"
        mv "$OUTPUT_DIR/${FILENAME}.fa.collapsed.map.filtered" "$OUTPUT_DIR/${FILENAME}.fa.collapsed.map"
    fi

    log "文件 $FILENAME 处理完成!"
    
    # 返回成功
    return 0
}

################################################################################
# 脚本参数设置及初始化
################################################################################

SCRIPT_PATH="$(cd "$(dirname "$0")" && pwd)" # 获取脚本所在目录的绝对路径
PROJECT_ROOT="$(dirname "$SCRIPT_PATH")" # 脚本在scripts目录下，所以项目根目录是脚本目录的父目录
LOG_DIR="$PROJECT_ROOT/logs"  # 默认日志目录
RESULTS_DIR="$PROJECT_ROOT/results"  # 默认结果目录

# 检查参数
if [ $# -lt 4 ]; then
    echo "用法: $0 <输入目录> <输出目录> <参考基因组路径> [--skip-existing] [--reverse]"
    exit 1
fi

# 获取必需的参数
INPUT_DIR="$1"
OUTPUT_DIR="$2"
PIRNA_REF="$3"

shift 3  # 移除已处理的参数

# 检查是否有--skip-existing参数和--reverse参数
SKIP_EXISTING=false
REVERSE_SEQ=false

# 处理剩余的参数
while [ $# -gt 0 ]; do
    if [ "$1" = "--skip-existing" ]; then
        SKIP_EXISTING=true
    elif [ "$1" = "--reverse" ]; then
        REVERSE_SEQ=true
    fi
    shift
done

################################################################################
# 目录和文件设置
################################################################################


TASK_ID=$(date +%Y%m%d_%H%M%S) # 获取当前时间作为任务ID
mkdir -p "$OUTPUT_DIR" # 创建输出目录、日志目录和结果目录
LOG_SUBDIR="$LOG_DIR/process_pirna_$TASK_ID" 
mkdir -p "$LOG_SUBDIR"  # 创建专门的日志子目录
RESULTS_SUBDIR="$RESULTS_DIR/pirna_analysis_$TASK_ID" 
mkdir -p "$RESULTS_SUBDIR" # 创建结果目录
MAIN_LOG="$LOG_SUBDIR/process_pirna.log"
touch "$MAIN_LOG" # 创建主日志文件


# 设置PYTHONPATH环境变量，确保能找到pirna_analysis模块
export PYTHONPATH="$PROJECT_ROOT/scripts:$PYTHONPATH"

log "开始处理piRNA数据"
log "脚本路径: $SCRIPT_PATH"
log "项目根目录: $PROJECT_ROOT"
log "输入目录: $INPUT_DIR"
log "输出目录: $OUTPUT_DIR"
log "日志目录: $LOG_SUBDIR"
log "结果目录: $RESULTS_SUBDIR"
log "参考基因组: $PIRNA_REF"
log "跳过已存在结果: $SKIP_EXISTING"
log "序列反向处理: $REVERSE_SEQ"

# 检查必要的脚本是否存在
SCRIPTS_DIR="$SCRIPT_PATH"
TBR2_DIR="$SCRIPTS_DIR/TBr2"
if [ ! -d "$TBR2_DIR" ]; then
    log "错误: TBr2脚本目录不存在: $TBR2_DIR"
    exit 1
fi

# 检查piRNA参考基因组是否存在
if [ ! -f "$PIRNA_REF" ]; then
    log "错误: piRNA参考基因组文件不存在: $PIRNA_REF"
    exit 1
fi

# 检查是否可以访问seqkit
if ! command -v seqkit &> /dev/null; then
    log "错误: seqkit未安装或未添加到PATH"
    exit 1
fi

# 检查是否可以访问seqmap
if ! command -v seqmap &> /dev/null; then
    log "错误: seqmap未安装或未添加到PATH"
    exit 1
fi

# 检查是否可以访问fastx_collapser
if ! command -v fastx_collapser &> /dev/null; then
    log "错误: fastx_collapser未安装或未添加到PATH"
    exit 1
fi

################################################################################
# 参考基因组分析
################################################################################

# 对piRNA参考基因组进行基本分析（只需要做一次）
REF_STATS_FILE="$OUTPUT_DIR/$(basename "$PIRNA_REF").stats"
if $SKIP_EXISTING && [ -f "$REF_STATS_FILE" ]; then
    log "跳过piRNA参考基因组分析: 结果文件已存在"
else
    log "对piRNA参考基因组进行基本分析..."
    perl "$TBR2_DIR/TBr2_basic-analyses.pl" -i "$PIRNA_REF" -o "$REF_STATS_FILE" >> "$MAIN_LOG" 2>&1
fi

# 如果需要反向处理，则生成反向参考基因组（只需要做一次）
if $REVERSE_SEQ; then
    # 将反向参考基因组保存在与原始参考基因组相同的文件夹中
    REF_DIR=$(dirname "$PIRNA_REF")
    REVERSE_REF="$REF_DIR/$(basename "$PIRNA_REF").reverse"
    
    if [ -f "$REVERSE_REF" ]; then
        log "跳过反向参考基因组生成: 文件已存在 ($REVERSE_REF)"
    else
        log "生成反向参考基因组: $REVERSE_REF"
        log "使用 20 个CPU核心生成反向参考基因组"
        seqkit seq -r "$PIRNA_REF" -j 20 -o "$REVERSE_REF" >> "$MAIN_LOG" 2>&1
        
        if [ $? -eq 0 ]; then
            log "反向参考基因组生成成功"
        else
            log "错误: 反向参考基因组生成失败"
            exit 1
        fi
    fi
fi

################################################################################
# 文件处理
################################################################################

# 查找所有fastq文件
FASTQ_FILES=()
for ext in "fq.gz" "fastq.gz"; do
    if ls "$INPUT_DIR"/*.$ext 1> /dev/null 2>&1; then
        for file in "$INPUT_DIR"/*.$ext; do
            FASTQ_FILES+=("$file")
        done
    fi
done

# 检查是否找到文件
if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
    log "错误: 在输入目录中没有找到fastq文件"
    exit 1
fi

log "找到 ${#FASTQ_FILES[@]} 个文件需要处理"

# 并行处理所有文件
PIDS=()
for FASTQ_FILE in "${FASTQ_FILES[@]}"; do
    # 启动新的处理进程，传递每个任务的CPU核心数
    process_file "$FASTQ_FILE" &
    PIDS+=($!)
    log "启动新进程 (PID: $!) 处理文件: $(basename "$FASTQ_FILE")"
done

# 等待所有进程完成
log "等待所有进程完成..."
for pid in "${PIDS[@]}"; do
    wait $pid
    log "进程 $pid 已完成"
done

log "所有文件处理完成!"

################################################################################
# 汇总分析piRNA分布
################################################################################

log "开始汇总分析piRNA分布..."
MAP_FILES=()
for FASTQ_FILE in "${FASTQ_FILES[@]}"; do
    # 获取原始文件名
    RAW_FILENAME=$(basename "$FASTQ_FILE" | sed 's/\.[^.]*$//')
    # 获取简化后的文件名
    SIMPLE_FILENAME=$(simplify_filename "$RAW_FILENAME")
    
    # 使用简化后的文件名查找map文件
    MAP_FILE="$OUTPUT_DIR/${SIMPLE_FILENAME}.fa.collapsed.map"
    if [ -f "$MAP_FILE" ]; then
        MAP_FILES+=("$MAP_FILE")
    else
        log "警告: 找不到map文件: $MAP_FILE (原始文件名: $RAW_FILENAME)"
    fi
done

if [ ${#MAP_FILES[@]} -gt 0 ]; then
    log "找到 ${#MAP_FILES[@]} 个map文件进行汇总分析"
    
    # 检查汇总分析结果是否已存在
    COMBINED_RESULT_EXISTS=false
    if $SKIP_EXISTING; then
        # 检查是否存在组合分布图（使用新的简化文件名格式）
        COMBINED_PLOT="$RESULTS_SUBDIR/distributions.png"
        COMBINED_UNIQUE_PLOT="$RESULTS_SUBDIR/unique_length_distribution.png"
        
        if [ -f "$COMBINED_PLOT" ] && [ -f "$COMBINED_UNIQUE_PLOT" ]; then
            COMBINED_RESULT_EXISTS=true
            log "跳过汇总分析: 结果文件已存在"
        fi
    fi
    
    if ! $COMBINED_RESULT_EXISTS; then
        # 构建命令行参数
        MAP_ARGS=""
        for MAP_FILE in "${MAP_FILES[@]}"; do
            # 将所有文件路径添加到一个-i参数后面
            if [ -z "$MAP_ARGS" ]; then
                MAP_ARGS="-i $MAP_FILE"
            else
                MAP_ARGS="$MAP_ARGS $MAP_FILE"
            fi
        done
        
        # 按照文件数量给analyze_pirna_distribution.py分配CPU核心数
        # 根据文件数量和可用核心数智能分配
        MAP_FILE_COUNT=${#MAP_FILES[@]}
        THREAD_ARGS="--threads $MAP_FILE_COUNT"
        ANALYSIS_CMD="python3 -m pirna_analysis.main $MAP_ARGS -o \"$RESULTS_SUBDIR\" -r $PIRNA_REF -v $THREAD_ARGS"
        log "执行命令: $ANALYSIS_CMD"
        
        # 执行汇总分析
        eval $ANALYSIS_CMD >> "$MAIN_LOG" 2>&1
    fi
    
    if [ $? -eq 0 ]; then
        log "piRNA分布汇总分析完成"
    else
        log "错误: piRNA分布汇总分析失败"
    fi
else
    log "错误: 没有找到可用的map文件进行汇总分析"
fi

log "详细日志文件保存在: $LOG_SUBDIR"

# 生成处理摘要
SUMMARY_FILE="$LOG_SUBDIR/summary.txt"
{
    echo "========== piRNA 处理摘要 ==========="
    echo "处理时间: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "脚本路径: $SCRIPT_PATH"
    echo "项目根目录: $PROJECT_ROOT"
    echo "输入目录: $INPUT_DIR"
    echo "输出目录: $OUTPUT_DIR"
    echo "日志目录: $LOG_SUBDIR"
    echo "结果目录: $RESULTS_SUBDIR"
    echo "参考基因组: $PIRNA_REF"
    echo "处理文件数: ${#FASTQ_FILES[@]}"
    echo "每个任务分配的CPU核心数: 10"
    echo "跳过已存在结果: $SKIP_EXISTING"
    echo "序列反向处理: $REVERSE_SEQ"
    echo "==============================="
    echo ""
    echo "处理的文件:"
    for file in "${FASTQ_FILES[@]}"; do
        echo "- $(basename "$file")"
    done
} > "$SUMMARY_FILE"

log "处理摘要已保存到: $SUMMARY_FILE"

# 创建结果目录的简单说明文件
README_FILE="$RESULTS_SUBDIR/README.txt"
{
    echo "========== piRNA 分析结果说明 ==========="
    echo "分析时间: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "输入目录: $INPUT_DIR"
    echo "处理文件数: ${#FASTQ_FILES[@]}"
    echo "日志目录: $LOG_SUBDIR"
    echo "中间文件目录: $OUTPUT_DIR"
    echo "==============================="
    echo ""
    echo "文件说明:"
    echo "- combined_pirna_distribution_*.png: 汇总的piRNA分布图"
    echo "- combined_pirna_distribution_*.tsv: 汇总的piRNA分布数据"
    echo ""
    echo "处理的文件:"
    for file in "${FASTQ_FILES[@]}"; do
        echo "- $(basename "$file")"
    done
} > "$README_FILE"

log "结果说明文件已保存到: $README_FILE"
