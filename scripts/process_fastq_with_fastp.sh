#!/bin/bash

# 脚本用于批量处理 fastq.gz/fq.gz 文件，使用 fastp 进行质量控制和过滤
# 作者：Cascade
# 创建日期：2025-03-24
# 更新日期：2025-05-13

# 使用方法
usage() {
    echo "用法: $0 [选项] 输入目录 输出目录"
    echo "选项:"
    echo "  -p | --paired       双端测序模式（默认为单端模式）"
    echo "  -h | --help         显示帮助信息"
    exit 1
}

# 处理命令行参数
PAIRED_MODE=false

while [ $# -gt 0 ]; do
    case "$1" in
        -p|--paired)
            PAIRED_MODE=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            if [ -z "$INPUT_DIR" ]; then
                INPUT_DIR="$1"
            elif [ -z "$OUTPUT_DIR" ]; then
                OUTPUT_DIR="$1"
            else
                echo "错误: 未知参数 $1"
                usage
            fi
            shift
            ;;
    esac
done

# 检查必要参数
if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "错误: 输入目录和输出目录是必需的"
    usage
fi

REPORT_DIR="$OUTPUT_DIR/fastp_reports"

# 创建输出目录（如果不存在）
mkdir -p "$OUTPUT_DIR"
mkdir -p "$REPORT_DIR"

# 输出处理信息
if [ "$PAIRED_MODE" = true ]; then
    echo "开始处理双端测序文件..."
else
    echo "开始处理单端测序文件（smRNA-seq模式）..."
fi

echo "输入目录: $INPUT_DIR"
echo "输出目录: $OUTPUT_DIR"
echo "报告目录: $REPORT_DIR"

# 计数器
if [ "$PAIRED_MODE" = true ]; then
    # 双端模式下计算R1文件数量
    total_files=$(find "$INPUT_DIR" \( -name "*_R1*.fastq.gz" -o -name "*_R1*.fq.gz" -o -name "*_1.fastq.gz" -o -name "*_1.fq.gz" \) | wc -l)
else
    # 单端模式下计算所有文件数量
    total_files=$(find "$INPUT_DIR" -name "*.fastq.gz" -o -name "*.fq.gz" | wc -l)
fi
processed=0

# 检查是否有文件需要处理
if [ $total_files -eq 0 ]; then
    echo "没有找到 fastq.gz 或 fq.gz 文件，退出处理。"
    exit 1
fi

# 根据模式处理文件
if [ "$PAIRED_MODE" = true ]; then
    # 处理双端测序文件
    # 找到所有R1文件
    for r1_file in $(find "$INPUT_DIR" \( -name "*_R1*.fastq.gz" -o -name "*_R1*.fq.gz" -o -name "*_1.fastq.gz" -o -name "*_1.fq.gz" \)); do
        if [ ! -f "$r1_file" ]; then
            continue
        fi

        # 确定R2文件命名方式
        if [[ "$r1_file" == *_R1*.fastq.gz ]]; then
            r2_file="${r1_file/_R1/_R2}"
        elif [[ "$r1_file" == *_R1*.fq.gz ]]; then
            r2_file="${r1_file/_R1/_R2}"
        elif [[ "$r1_file" == *_1.fastq.gz ]]; then
            r2_file="${r1_file/_1/_2}"
        elif [[ "$r1_file" == *_1.fq.gz ]]; then
            r2_file="${r1_file/_1/_2}"
        fi

        # 检查R2文件是否存在
        if [ ! -f "$r2_file" ]; then
            echo "警告: 找不到对应的R2文件: $r2_file, 跳过处理 $r1_file"
            continue
        fi

        # 获取文件名（不含路径和扩展名）
        if [[ "$r1_file" == *.fastq.gz ]]; then
            base_name=$(basename "$r1_file" | sed 's/_R1\(.*\)\.fastq\.gz$\|_1\(.*\)\.fastq\.gz$//')
        else
            base_name=$(basename "$r1_file" | sed 's/_R1\(.*\)\.fq\.gz$\|_1\(.*\)\.fq\.gz$//')
        fi

        echo "处理双端文件: $base_name (第 $((processed+1)) / $total_files 对)"

        # 输出文件路径
        output_r1="$OUTPUT_DIR/${base_name}_R1.clean.fastq.gz"
        output_r2="$OUTPUT_DIR/${base_name}_R2.clean.fastq.gz"

        # fastp 报告文件
        html_report="$REPORT_DIR/${base_name}.html"
        json_report="$REPORT_DIR/${base_name}.json"

        # 运行 fastp - 双端测序参数
        fastp -i "$r1_file" -I "$r2_file" -o "$output_r1" -O "$output_r2" \
            --detect_adapter_for_pe --complexity_threshold 30 --thread 16 \
            --cut_mean_quality 20 --qualified_quality_phred 15 --unqualified_percent_limit 40 \
            --n_base_limit 5 --length_required 50 --low_complexity_filter \
            --html "$html_report" --json "$json_report"

        # 检查 fastp 是否成功
        if [ $? -eq 0 ]; then
            processed=$((processed + 1))
            echo "成功处理: $base_name"
        else
            echo "处理失败: $base_name"
        fi

        echo "进度: $processed / $total_files"
        echo "-----------------------------------"
    done
else
    # 处理单端测序文件
    for fastq_file in "$INPUT_DIR"/*.fastq.gz "$INPUT_DIR"/*.fq.gz; do
        # 检查文件是否存在和是否是常规文件
        if [ ! -f "$fastq_file" ]; then
            continue
        fi

        # 获取文件名（不含路径和扩展名）
        if [[ "$fastq_file" == *.fastq.gz ]]; then
            filename=$(basename "$fastq_file" .fastq.gz)
        else
            filename=$(basename "$fastq_file" .fq.gz)
        fi

        # 输出文件路径
        output_file="$OUTPUT_DIR/${filename}.clean.fastq.gz"

        # fastp 报告文件
        html_report="$REPORT_DIR/${filename}.html"
        json_report="$REPORT_DIR/${filename}.json"

        echo "处理文件: $filename"

        # 运行 fastp - smRNA-seq单端测序参数
        # Step 1: Remove adapter and perform initial filtering
        fastp -i "$fastq_file" -o "${OUTPUT_DIR}/${filename}.tmp1.fq.gz" \
            -a AACTGTAGGCACCATCAAT --length_required 1 --cut_front -W 4 -M 1 -q 1 -n 20 -w 16

        # Step 2: Perform additional filtering without adapter trimming
        fastp -i "${OUTPUT_DIR}/${filename}.tmp1.fq.gz" -o "${OUTPUT_DIR}/${filename}.tmp2.fq.gz" \
            -A -t 0 -w 16 --length_required 1 --cut_front -W 4 -M 1 -q 1 -n 20

        # Step 3: Final filtering and generate reports
        fastp -i "${OUTPUT_DIR}/${filename}.tmp2.fq.gz" -o "$output_file" \
            -A --length_required 18 --length_limit 37 -w 16 -n 0 -z 9 -q 20 -W 4 -M 20 \
            --html "$html_report" --json "$json_report"

        # Clean up temporary files
        rm "${OUTPUT_DIR}/${filename}.tmp1.fq.gz" "${OUTPUT_DIR}/${filename}.tmp2.fq.gz"


        # 检查 fastp 是否成功
        if [ $? -eq 0 ]; then
            processed=$((processed + 1))
            echo "成功处理: $filename"
        else
            echo "处理失败: $filename"
        fi

        echo "进度: $processed / $total_files"
        echo "-----------------------------------"
    done
fi

if [ "$PAIRED_MODE" = true ]; then
    echo "处理完成！共处理 $processed 对双端测序文件。"
else
    echo "处理完成！共处理 $processed 个单端测序文件。"
fi
