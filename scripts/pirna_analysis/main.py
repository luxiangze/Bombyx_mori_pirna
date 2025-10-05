#!/usr/bin/env python3
"""
主程序模块 - 负责命令行参数解析和程序流程控制
"""

import os
import sys
import time
import argparse
import multiprocessing as mp
from tqdm import tqdm
import pandas as pd

# 导入自定义模块
from pirna_analysis.data_parser import parse_map_file_optimized
from pirna_analysis.data_analyzer import (
    analyze_length_distribution,
    analyze_5prime_shift,
    analyze_3prime_position,
    analyze_unique_reads_length,
    analyze_base_composition_by_position
)
from pirna_analysis.visualization import (
    plot_combined_distributions,
    plot_combined_unique_reads_length,
    plot_base_composition_logo
)
from pirna_analysis.utils import extract_sample_name

def process_file(map_file, ref_file, output_prefix, verbose=False):
    """
    处理单个map文件的函数，用于多进程并行处理
    只负责数据处理和保存，不进行绘图
    """
    sample_name = extract_sample_name(map_file)
    if verbose:
        print(f"处理样本 {sample_name}: {map_file}")
    
    # 使用优化版本的解析函数
    read_info, read_to_pirna_count, read_to_pirnas, filter_stats = parse_map_file_optimized(map_file, ref_file)
    
    # 碱基组成分析
    df_base = analyze_base_composition_by_position(read_info, read_to_pirna_count, read_to_pirnas)
    
    if verbose:
        total_reads, filtered_reads = filter_stats
        print(f"  在{sample_name}中找到{len(read_info)}个唯一reads")
        print(f"  过滤了{filtered_reads}/{total_reads}个起始位置>9的reads（调整后>6）")
    
    # 分析长度分布
    df_length = analyze_length_distribution(read_info, read_to_pirna_count)
    
    # 分析5'端位移分布
    df_5prime = analyze_5prime_shift(read_info, read_to_pirna_count, read_to_pirnas)
    
    # 分析3'端相对位置分布
    df_3prime = analyze_3prime_position(read_info, read_to_pirna_count, read_to_pirnas)
    
    # 分析去重后的reads长度分布
    df_unique = analyze_unique_reads_length(read_info)
    
    # 创建输出目录
    os.makedirs(output_prefix, exist_ok=True)
    
    # 保存结果
    df_length.to_csv(f"{output_prefix}/{sample_name}_length_distribution.csv")
    df_5prime.to_csv(f"{output_prefix}/{sample_name}_5prime_shift.csv")
    df_3prime.to_csv(f"{output_prefix}/{sample_name}_3prime_position.csv")
    df_unique.to_csv(f"{output_prefix}/{sample_name}_unique_reads_length.csv")
    df_base.to_csv(f"{output_prefix}/{sample_name}_base_composition.csv")
    
    # 返回样本名称和处理结果
    return sample_name

def main():
    # 设置命令行参数解析
    parser = argparse.ArgumentParser(description='分析piRNA的长度分布、5\'端位移分布、3\'端相对位置分布和去重后的reads长度分布')
    parser.add_argument('-i', '--input', nargs='+', required=True, help='输入的map文件路径，可以指定多个文件')
    parser.add_argument('-o', '--output', required=True, help='输出目录或前缀')
    parser.add_argument('-v', '--verbose', action='store_true', help='显示详细处理信息')
    parser.add_argument('-p', '--prefix', default='', help='输出文件的前缀，默认为空')
    parser.add_argument('-r', '--ref', default='/home/gyk/project/ld_pirna/data/bmo.v3.0.fa', help='参考piRNA文件路径')
    parser.add_argument('--threads', type=int, default=0, help='并行处理的线程数，默认为0表示使用所有可用CPU核心')
    
    # 解析命令行参数
    args = parser.parse_args()
    
    # 获取输入文件列表和输出前缀
    map_files = args.input
    output_prefix = args.output
    verbose = args.verbose
    file_prefix = args.prefix
    
    # 设置多进程的线程数
    threads = args.threads
    if threads <= 0:
        threads = mp.cpu_count()
    
    # 如果输出目录不存在，则创建
    output_dir = os.path.dirname(output_prefix) if os.path.dirname(output_prefix) else output_prefix
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    print(f"使用 {threads} 个CPU核心进行处理")
    
    # 处理多个样本
    all_data = {}
    all_unique_data = {}
    start_time = time.time()
    
    # 使用多进程并行处理文件
    if len(map_files) > 1 and threads > 1:
        print(f"并行处理 {len(map_files)} 个文件...")
        # 创建进程池
        pool = mp.Pool(min(threads, len(map_files)))
        # 准备处理函数的参数
        process_args = [(map_file, args.ref, output_prefix, verbose) for map_file in map_files]
        sample_names = pool.starmap(process_file, process_args)
        pool.close()
        pool.join()
    else:
        sample_names = []
        for map_file in map_files:
            sample_name = process_file(map_file, args.ref, output_prefix, verbose)
            sample_names.append(sample_name)

    # 汇总结果，读取各类csv
    for sample_name in sample_names:
        sample_prefix = os.path.join(args.output, sample_name)
        df_length = pd.read_csv(f"{sample_prefix}_length_distribution.csv")
        df_5prime = pd.read_csv(f"{sample_prefix}_5prime_shift.csv")
        df_3prime = pd.read_csv(f"{sample_prefix}_3prime_position.csv")
        df_unique = pd.read_csv(f"{sample_prefix}_unique_reads_length.csv")
        df_base = pd.read_csv(f"{sample_prefix}_base_composition.csv", index_col=0)
        all_data[sample_name] = (df_length, df_5prime, df_3prime)
        all_unique_data[sample_name] = df_unique
        # 绘制logo
        plot_base_composition_logo(df_base, f"{sample_prefix}_base_logo.pdf", title=f"{sample_name} base composition")

    # 始终绘制组合分布图，即使只有一个样本
    if verbose:
        print("生成组合分布图...")
    
    # 绘制常规分布的组合图，将三个子图分别输出为单独的图片
    plot_combined_distributions(all_data, f"{output_prefix}/{file_prefix}")
    
    # 绘制去重后reads长度分布的组合图
    plot_combined_unique_reads_length(all_unique_data, f"{output_prefix}/{file_prefix}")
    
    print(f"分析完成。结果已保存到: {output_prefix}")
    print(f"总耗时: {time.time() - start_time:.2f} 秒")

if __name__ == "__main__":
    main()
