# piRNA Analysis Scripts Guide

This directory contains various scripts for analyzing piRNA data. These scripts are used to process sequencing data, analyze piRNA length distributions and positional features, and identify unprocessed sequences.

## Quick Start

Below is a complete example workflow for piRNA analysis:

```bash
# 1. Process raw sequencing data
./process_pirna.sh /path/to/fastq_files /path/to/output 8

# 2. Analyze piRNA length distribution and positional features
./batch_analyze_pirna.sh /path/to/output/map_files /path/to/results Control.map

# 3. Plot sample-to-sample comparison using the R script
Rscript plot_scatter_comparison.R --input_dir=/path/to/results --output_dir=/path/to/figures \
  --sample_config=/path/to/sample_config.txt
```

If you only need to replot figures from existing CSV results, use:

```bash
# Replot figures from previously generated CSV files
python3 replot_from_csv.py -i /path/to/csv_files -o /path/to/replotted_figures

# Or plot sample-to-sample comparison using the R script
Rscript plot_scatter_comparison.R --input_dir=/path/to/csv_files --output_dir=/path/to/figures
```

## 主要脚本

### 数据处理脚本

#### `process_pirna.sh`

**功能**：自动处理 piRNA 数据的主脚本，包括数据预处理、比对和分析。

**用法**：
```bash
./process_pirna.sh <输入目录> <输出目录> [并行任务数] [参考基因组路径] [日志目录] [--skip-existing]
```

**参数说明**：
- `<输入目录>`：包含原始 fastq 文件的目录
- `<输出目录>`：处理结果的输出目录
- `[并行任务数]`：可选，同时处理的文件数量，默认为 4
- `[参考基因组路径]`：可选，piRNA 参考基因组文件路径
- `[日志目录]`：可选，日志文件保存目录
- `[--skip-existing]`：可选，如果结果文件已存在则跳过处理

#### `process_fastq_with_fastp.sh`

**功能**：使用 fastp 工具处理 fastq 文件，进行质量控制和过滤。

**用法**：
```bash
./process_fastq_with_fastp.sh <输入文件> <输出目录>
```

### piRNA 分析模块

#### `pirna_analysis` 包

这是一个模块化的 Python 包，用于分析 piRNA 的长度分布、5'端位移分布和3'端相对位置分布。包含以下模块：

- `data_parser.py`：负责解析 map 文件和参考 piRNA 文件
- `data_analyzer.py`：负责数据分析功能
- `visualization.py`：负责绘图功能
- `utils.py`：提供辅助功能
- `main.py`：主程序入口和流程控制

#### `analyze_pirna.py`

**功能**：piRNA 分析的入口脚本，调用 pirna_analysis 包进行分析。

**用法**：
```bash
python3 analyze_pirna.py -i <map_file1> [<map_file2> ...] -o <output_directory> [其他参数]
```

**参数说明**：
- `-i, --input`：输入的 map 文件路径，可以指定多个文件
- `-o, --output`：输出目录
- `-v, --verbose`：显示详细处理信息
- `-p, --prefix`：输出文件的前缀
- `-r, --ref`：参考 piRNA 文件路径
- `--threads`：并行处理的线程数

#### `replot_from_csv.py`

**功能**：从已生成的CSV文件重新绘制piRNA分析图表，避免重新解析map文件。

**特点**：
1. 自动查找指定目录中的所有CSV结果文件
2. 支持绘制四种图表：长度分布、5'端位移分布、3'端相对位置分布和唯一读数长度分布
3. 确保所有图表的横坐标都显示为整数
4. 可以指定输出目录和文件前缀

**用法**：
```bash
python3 replot_from_csv.py -i <CSV文件目录> -o <输出目录> [-p <输出文件前缀>]
```

**参数说明**：
- `-i, --input`：包含CSV文件的输入目录
- `-o, --output`：输出目录
- `-p, --prefix`：可选，输出文件前缀

#### `analyze_pirna_distribution.py`（已拆分为模块）

**功能**：分析 piRNA 的长度分布、5'端位移分布和3'端相对位置分布。

**注意**：此脚本已被拆分为 `pirna_analysis` 包，建议使用 `analyze_pirna.py` 或直接调用 `pirna_analysis.main` 模块。

### 批处理脚本

#### `batch_analyze_pirna.sh`

**功能**：批量分析指定目录中的所有 map 文件，自动分配 CPU 资源。

**特点**：
1. 自动查找指定目录中的所有 .map 文件
2. 允许指定控制组文件（默认为 Control.fa.collapsed.no-dust.map）
3. 自动计算并合理分配 CPU 资源（使用服务器 90% 的 CPU）
4. 同时分析 32nt 和 28nt 的未处理序列
5. 生成带时间戳的日志文件
6. 提供详细的运行信息和错误处理

### 序列分析脚本

#### `identify_unprocessed_sequences.py`

**功能**：识别未从 32nt 处理到 28nt 的序列以及未从 28nt 处理到其他长度的序列。

**特点**：
1. 专注于一对样本的分析
2. 允许用户自定义样本名称
3. 优化的输出结构

**用法**：
```bash
python3 identify_unprocessed_sequences.py <对照组文件> <实验组文件> -o <输出目录> [其他参数]
```

#### `plot_length_distribution.py`

**功能**：绘制 piRNA 长度分布图。

### 辅助脚本

#### `count_length.sh`

**功能**：统计序列长度分布。

#### `rename.sh`

**功能**：批量重命名文件。

#### `piPipes_run.sh` 和 `tailor_run.sh`

**功能**：调用 piPipes 和 Tailor 工具进行 piRNA 分析的脚本。

## 使用示例

### 基本处理流程

1. 使用 `process_pirna.sh` 处理原始数据：
   ```bash
   ./process_pirna.sh /path/to/fastq_files /path/to/output 8
   ```

2. 使用 `analyze_pirna.py` 分析处理后的数据：
   ```bash
   python3 analyze_pirna.py -i /path/to/output/*.map -o /path/to/results -v
   ```

3. 使用 `identify_unprocessed_sequences.py` 识别未处理序列：
   ```bash
   python3 identify_unprocessed_sequences.py /path/to/control.map /path/to/treatment.map -o /path/to/results
   ```

4. 使用 `replot_from_csv.py` 从已生成的CSV文件重新绘制图表：
   ```bash
   python3 replot_from_csv.py -i /path/to/csv_files -o /path/to/replotted_figures -p "replotted_"
   ```

5. 使用 `plot_scatter_comparison.R` 绘制样本间的散点图比较：
   ```bash
   Rscript plot_scatter_comparison.R --input_dir=/path/to/csv_files --output_dir=/path/to/output
   ```

### 批量处理

使用 `batch_analyze_pirna.sh` 批量分析多个样本：
```bash
./batch_analyze_pirna.sh /path/to/map_files /path/to/output Control.map
```

### 数据可视化

#### `plot_scatter_comparison.R`

**功能**：通用散点图比较工具，用于绘制两个样本之间的散点图比较，特别适合piRNA剪切数据的比较。使用ggplot2绘制高质量的散点图，支持多种自定义选项。

**特点**：
1. 自动检查并安装所需的R包
2. 自动处理多个样本，生成所有样本对之间的比较图
3. 动态调整颜色映射范围，根据数据的fold change分布
4. 支持对数坐标轴和线性坐标轴
5. 美观的可视化效果，包括对角线参考线和颜色条
6. 支持样本配置文件，自动识别对照组和处理组
7. 根据样本类型自动调整样本顺序，确保对照组在x轴，处理组在y轴

**用法**：
```bash
Rscript plot_scatter_comparison.R --input_dir=/path/to/csv_files --output_dir=/path/to/output [--sample_config=/path/to/sample_config.txt]
```

**参数说明**：
- `--input_dir`：输入CSV文件所在目录（必需）
- `--output_dir`：输出PDF文件保存目录（必需）
- `--sample_config`：样本配置文件路径，包含样本名称和类型信息（可选）
- `--pattern`：用于筛选CSV文件的模式，默认为"_pirna_trim.csv"
- `--id_column`：用于合并数据的ID列名，默认为"piRNA_ID"
- `--value_column`：用于比较的数值列名，默认为"Trim_index"
- `--log_scale`：是否使用对数坐标轴，默认为TRUE
- `--point_size`：点的大小，默认为1
- `--point_alpha`：点的透明度，默认为0.7
- `--width`：输出PDF的宽度(英寸)，默认为7
- `--height`：输出PDF的高度(英寸)，默认为7

**样本配置文件格式**：
样本配置文件是一个简单的文本文件，每行包含一个样本的信息，格式为：
```
<样本名称> <样本类型>
```
其中样本类型可以是`untreatment`（对照组）或`treatment`（处理组）。

示例：
```
Control untreatment
DHX15-KD3 treatment
SUGP1-KD2 treatment
RNPS1-KD1 treatment
```

**示例**：
```bash
# 基本用法
Rscript plot_scatter_comparison.R --input_dir=/path/to/csv_files --output_dir=/path/to/output

# 使用样本配置文件
Rscript plot_scatter_comparison.R --input_dir=/path/to/csv_files --output_dir=/path/to/output \
  --sample_config=/path/to/sample_config.txt

# 自定义参数
Rscript plot_scatter_comparison.R --input_dir=/path/to/csv_files --output_dir=/path/to/output \
  --point_size=1.5 --point_alpha=0.5 --log_scale=FALSE --width=10 --height=8

# 使用不同的文件模式和列名
Rscript plot_scatter_comparison.R --input_dir=/path/to/csv_files --output_dir=/path/to/output \
  --pattern="_expression.csv" --id_column="gene_id" --value_column="FPKM"
```

#### `plot_featrue_distribution.R`

**功能**：使用tidyplots包绘制环状图（Donut Plot），用于可视化不同类别的分布情况。

**特点**：
1. 基于tidyplots包实现，生成美观的环状图
2. 支持自定义颜色方案、图表风格、标签大小等
3. 自动计算并显示每个类别的百分比
4. 支持浅色和深色两种主题
5. 可调整环形宽度和类别排序

**用法**：
```bash
Rscript plot_featrue_distribution.R -i <CSV文件路径> -o <输出图片路径> [其他参数]
```

**参数说明**：
- `-i, --input`：输入CSV文件路径，第一列为类别，第二列为值（必需）
- `-o, --output`：输出图片文件路径，支持PDF、PNG等格式（默认：donut_plot.pdf）
- `-t, --title`：图表标题（默认："特征分布"）
- `-w, --width`：输出图片宽度（英寸）（默认：8）
- `-e, --height`：输出图片高度（英寸）（默认：8）
- `-c, --color_palette`：颜色调色板名称，选项有friendly, neutral, vibrant等（默认：friendly）
- `-l, --label_size`：标签文字大小（默认：12）
- `-p, --percentage`：是否显示百分比（默认：TRUE）
- `-d, --donut_width`：环状图的宽度，范围为0-1（默认：0.6）
- `-s, --style`：图表风格，选项有light或dark（默认：light）
- `-r, --reverse`：是否反转类别顺序（默认：FALSE）

**示例**：
```bash
# 基本用法
Rscript plot_featrue_distribution.R -i data/feature_counts.csv -o figures/feature_dist.pdf

# 自定义参数
Rscript plot_featrue_distribution.R -i data/feature_counts.csv -o figures/feature_dist.pdf \
  -t "基因特征分布" -c vibrant -s dark -d 0.8 -l 14 -r
```

**CSV文件格式示例**：
```
类别,值
基因启动子,150
基因终止子,120
基因间区,200
外显子,300
内含子,250
```

## 注意事项

1. 大多数脚本依赖于特定的文件命名格式，请确保输入文件的命名符合要求。
2. 处理大量数据时，请注意合理分配 CPU 资源，避免服务器过载。
3. 使用 `--skip-existing` 参数可以跳过已处理的文件，节省时间。
4. 使用R脚本（如`plot_scatter_comparison.R`）需要安装R环境，脚本会自动安装所需的R包。
5. 对于piRNA剪切数据的比较，建议使用R脚本`plot_scatter_comparison.R`而非Python脚本，因为R脚本能够正确地显示所有数据点。
6. 当处理大量数据点（如每个样本40万+数据点）时，生成的PDF文件可能会较大（约8-10MB），请确保有足够的磁盘空间。

#### process_fastq_with_fastp.sh

```bash
# 单端模式（适用于smRNA-seq）
scripts/process_fastq_with_fastp.sh 输入目录 输出目录

# 双端模式（适用于常规双端测序）
scripts/process_fastq_with_fastp.sh -p 输入目录 输出目录

# 查看帮助
scripts/process_fastq_with_fastp.sh -h
```