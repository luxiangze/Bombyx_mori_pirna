#!/bin/bash

# 批量分析piRNA序列处理情况的脚本
# 用法: ./batch_analyze_pirna.sh <map文件目录> <输出目录> <配置文件>

# 检查参数
if [ $# -lt 3 ]; then
    echo "用法: $0 <map文件目录> <输出目录> <配置文件>"
    echo "示例: $0 /path/to/maps /path/to/output /path/to/config.txt"
    echo ""
    echo "配置文件格式示例:"
    echo "Control untreatment"
    echo "DHX15-KD3 treatment"
    echo "SUGP1-KD2 treatment"
    echo "RNPS1-KD1 treatment"
    exit 1
fi

# 获取参数
MAP_DIR="$1"
OUTPUT_DIR="$2"
CONFIG_FILE="$3"
# 如果没有指定日志目录，则使用输出目录
LOG_DIR="${4:-$OUTPUT_DIR}"

# 检查目录和文件是否存在
if [ ! -d "$MAP_DIR" ]; then
    echo "错误: 目录 $MAP_DIR 不存在"
    exit 1
fi

if [ ! -f "$CONFIG_FILE" ]; then
    echo "错误: 配置文件 $CONFIG_FILE 不存在"
    exit 1
fi

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

# 读取配置文件，获取对照组和实验组
declare -A SAMPLE_TYPES
CONTROL_SAMPLES=()
TREATMENT_SAMPLES=()

while read -r sample type; do
    SAMPLE_TYPES["$sample"]="$type"
    if [ "$type" == "untreatment" ]; then
        CONTROL_SAMPLES+=("$sample")
    elif [ "$type" == "treatment" ]; then
        TREATMENT_SAMPLES+=("$sample")
    else
        echo "警告: 未知的样本类型 '$type'，样本 '$sample' 将被忽略"
    fi
done < "$CONFIG_FILE"

# 检查是否找到对照组和实验组
if [ ${#CONTROL_SAMPLES[@]} -eq 0 ]; then
    echo "错误: 配置文件中未找到对照组样本 (untreatment)"
    exit 1
fi

if [ ${#TREATMENT_SAMPLES[@]} -eq 0 ]; then
    echo "错误: 配置文件中未找到实验组样本 (treatment)"
    exit 1
fi

echo "从配置文件中读取到:"
echo "  对照组样本: ${CONTROL_SAMPLES[*]}"
echo "  实验组样本: ${TREATMENT_SAMPLES[*]}"

# 查找所有map文件
declare -A MAP_FILES

# 首先列出所有map文件
echo "在 $MAP_DIR 中找到的map文件:"
while IFS= read -r -d '' file; do
    echo "  $(basename "$file")"
done < <(find "$MAP_DIR" -name "*.map" -type f -print0)

# 然后匹配样本名称
while IFS= read -r -d '' file; do
    filename=$(basename "$file")
    for sample in "${!SAMPLE_TYPES[@]}"; do
        if [[ "$filename" == *"$sample"* ]]; then
            MAP_FILES["$sample"]="$file"
            echo "匹配: 样本 '$sample' -> 文件 '$(basename "$file")'" 
            break
        fi
    done
done < <(find "$MAP_DIR" -name "*.map" -type f -print0)

# 检查是否找到所有样本的map文件
MISSING_SAMPLES=()
for sample in "${!SAMPLE_TYPES[@]}"; do
    if [ -z "${MAP_FILES[$sample]}" ]; then
        MISSING_SAMPLES+=("$sample")
    fi
done

if [ ${#MISSING_SAMPLES[@]} -gt 0 ]; then
    echo "警告: 以下样本未找到对应的map文件:"
    for sample in "${MISSING_SAMPLES[@]}"; do
        echo "  $sample"
    done
    
    echo ""
    echo "可用的map文件:"
    while IFS= read -r -d '' file; do
        echo "  $(basename "$file")"
    done < <(find "$MAP_DIR" -name "*.map" -type f -print0)
    
    echo ""
    echo "是否继续分析? [y/N]"
    read -r answer
    if [[ ! "$answer" =~ ^[Yy]$ ]]; then
        echo "分析已取消"
        exit 1
    fi
fi

# 创建时间戳和日志文件
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
# 确保日志目录存在
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/analysis_log_$TIMESTAMP.log"

# 对每个对照组样本进行分析
for control_sample in "${CONTROL_SAMPLES[@]}"; do
    control_file="${MAP_FILES[$control_sample]}"
    
    # 如果没有找到对照组文件，跳过
    if [ -z "$control_file" ]; then
        echo "警告: 未找到对照组样本 $control_sample 的map文件，跳过"
        continue
    fi
    
    echo "使用 $control_sample 作为对照组进行分析..." | tee -a "$LOG_FILE"
    
    # 对每个实验组样本进行分析
    for treatment_sample in "${TREATMENT_SAMPLES[@]}"; do
        treatment_file="${MAP_FILES[$treatment_sample]}"
        
        # 如果没有找到实验组文件，跳过
        if [ -z "$treatment_file" ]; then
            echo "  警告: 未找到实验组样本 $treatment_sample 的map文件，跳过" | tee -a "$LOG_FILE"
            continue
        fi
        
        echo "  分析 $control_sample vs $treatment_sample..." | tee -a "$LOG_FILE"
        
        # 创建样本特定的输出目录，使用简洁的名称
        sample_output_dir="$OUTPUT_DIR/${control_sample}_vs_${treatment_sample}"
        mkdir -p "$sample_output_dir"
        
        # 运行32nt和28nt的分析
        cmd="python $(dirname "$0")/identify_unprocessed_sequences.py -i \"$control_file\" \"$treatment_file\" -o \"$sample_output_dir\" -t both -v --control-name \"$control_sample\" --exp-name \"$treatment_sample\""
        echo "  $cmd" | tee -a "$LOG_FILE"
        
        # 执行命令
        python "$(dirname "$0")/identify_unprocessed_sequences.py" -i "$control_file" "$treatment_file" -o "$sample_output_dir" -t both -v --control-name "$control_sample" --exp-name "$treatment_sample" 2>&1 | tee -a "$LOG_FILE"
        
        # 重命名输出目录中可能的冗长文件夹名称
        for dir in "$sample_output_dir"/*; do
            if [ -d "$dir" ]; then
                dir_name=$(basename "$dir")
                if [[ "$dir_name" == *"_vs_"* ]]; then
                    # 已经是简洁名称，不需要重命名
                    continue
                fi
                # 创建简洁的目录名
                new_dir_name="${control_sample}_vs_${treatment_sample}"
                if [ "$dir_name" != "$new_dir_name" ]; then
                    mv "$dir" "$sample_output_dir/$new_dir_name"
                fi
            fi
        done
        
        echo "  完成 $control_sample vs $treatment_sample 的分析" | tee -a "$LOG_FILE"
        echo "" | tee -a "$LOG_FILE"
    done
done

echo "所有分析完成，结果保存在 $OUTPUT_DIR"
echo "日志文件: $LOG_FILE"
