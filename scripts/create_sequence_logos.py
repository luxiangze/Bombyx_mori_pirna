#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
分析CSV文件中的序列特征并生成序列logo图

此脚本用于分析输入文件夹中以_length.csv为后缀的文件，
提取第一列的序列数据，使用logomaker生成序列logo图，
并将结果保存为PDF格式。
'''

import os
import sys
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logomaker
import argparse
from collections import Counter
import re

def parse_args():
    '''
    解析命令行参数
    '''
    parser = argparse.ArgumentParser(description='分析序列特征并生成序列logo图')
    parser.add_argument('-i', '--input_dir', required=True, help='输入文件夹路径，包含以_length.csv为后缀的文件')
    parser.add_argument('-l', '--length_filter', type=int, default=None, help='可选参数，只分析特定长度的序列')
    parser.add_argument('-n', '--top_n', type=int, default=100, help='用于生成logo的序列数量，默认为前100个序列')
    parser.add_argument('-p', '--position', type=str, default='all', choices=['all', '5p', '3p'], help='分析位置，可选all（全序列）、5p（5\'端）或3p（3\'端）')
    parser.add_argument('--logfc_filter', type=float, default=None, help='可选参数，根据logFC值筛选序列，例如">2"表示logFC>2的序列')
    return parser.parse_args()

def create_position_frequency_matrix(sequences, position='all'):
    '''
    创建位置频率矩阵
    
    参数:
    sequences - 序列列表
    position - 分析位置，'all'表示全序列，'5p'表示5'端，'3p'表示3'端
    
    返回:
    位置频率矩阵
    '''
    # 确定序列长度
    if position == 'all':
        # 使用所有序列
        pass
    elif position == '5p':
        # 只使用序列的前10个碱基（如果序列长度足够）
        sequences = [seq[:10] if len(seq) >= 10 else seq for seq in sequences]
    elif position == '3p':
        # 只使用序列的后10个碱基（如果序列长度足够）
        sequences = [seq[-10:] if len(seq) >= 10 else seq for seq in sequences]
    
    # 找到最长序列的长度
    max_length = max(len(seq) for seq in sequences)
    
    # 初始化计数矩阵
    counts_matrix = {}
    for base in 'ACGTU':
        counts_matrix[base] = [0] * max_length
    
    # 计算每个位置的碱基频率
    for seq in sequences:
        for i, base in enumerate(seq):
            base = base.upper()
            if base in counts_matrix:
                counts_matrix[base][i] += 1
    
    # 转换为DataFrame
    counts_df = pd.DataFrame(counts_matrix)
    
    # 将'U'合并到'T'中（如果存在）
    if 'U' in counts_df.columns:
        if 'T' in counts_df.columns:
            counts_df['T'] = counts_df['T'] + counts_df['U']
        else:
            counts_df['T'] = counts_df['U']
        counts_df = counts_df.drop('U', axis=1)
    
    # 计算每个位置的总计数
    totals = counts_df.sum(axis=1)
    
    # 将计数转换为频率
    freq_matrix = counts_df.div(totals, axis=0).fillna(0)
    
    return freq_matrix

def create_sequence_logo(freq_matrix, output_file, title=None):
    '''
    使用logomaker创建序列logo图
    
    参数:
    freq_matrix - 位置频率矩阵
    output_file - 输出文件路径
    title - 图表标题
    '''
    # 设置中文字体支持
    plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans', 'Bitstream Vera Sans', 'sans-serif']
    plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题
    
    plt.figure(figsize=(10, 3))
    
    # 创建Logo对象，使用系统默认字体
    logo = logomaker.Logo(freq_matrix, shade_below=.5, fade_below=.5)
    
    # 样式设置
    logo.style_spines(visible=False)
    logo.style_xticks(rotation=0, fmt='%d', anchor=0)
    
    # 设置y轴范围和标签
    plt.ylim(0, 1)
    plt.ylabel('Frequency')
    plt.xlabel('Position')
    
    # 设置标题
    if title:
        plt.title(title)
    
    # 保存图像
    plt.tight_layout()
    plt.savefig(output_file, format='pdf', dpi=300)
    plt.close()
    
    print(f"已生成序列logo图: {output_file}")

def process_csv_file(file_path, output_dir, length_filter=None, top_n=100, position='all', logfc_filter=None):
    '''
    处理CSV文件并生成序列logo图
    
    参数:
    file_path - CSV文件路径
    output_dir - 输出目录
    length_filter - 序列长度过滤条件
    top_n - 用于生成logo的序列数量
    position - 分析位置
    logfc_filter - logFC过滤条件
    '''
    print(f"\n处理文件: {os.path.basename(file_path)}")
    try:
        # 读取CSV文件
        df = pd.read_csv(file_path)
        
        # 检查文件格式
        if len(df.columns) < 1:
            print(f"警告: {file_path} 格式不正确，至少需要1列")
            return
        
        # 获取第一列（序列列）
        sequences = df.iloc[:, 0].tolist()
        
        # 如果存在logFC列，根据logFC值筛选序列
        if logfc_filter is not None and 'logFC' in df.columns:
            match = re.match(r'([<>])([\d.-]+)', str(logfc_filter))
            if match:
                operator, value = match.groups()
                value = float(value)
                if operator == '>':
                    df = df[df['logFC'] > value]
                elif operator == '<':
                    df = df[df['logFC'] < value]
                sequences = df.iloc[:, 0].tolist()
        
        # 如果存在Length列，根据长度筛选序列
        if length_filter is not None and 'Length' in df.columns:
            df = df[df['Length'] == length_filter]
            sequences = df.iloc[:, 0].tolist()
        
        # 确保我们有足够的序列
        if len(sequences) == 0:
            print(f"警告: {file_path} 在应用过滤条件后没有剩余序列")
            return
        
        # 只使用前N个序列
        sequences = sequences[:min(top_n, len(sequences))]
        
        # 创建位置频率矩阵
        freq_matrix = create_position_frequency_matrix(sequences, position)
        
        # 生成输出文件名
        base_name = os.path.basename(file_path).replace('_length.csv', '')
        position_suffix = '' if position == 'all' else f'_{position}'
        length_suffix = '' if length_filter is None else f'_length{length_filter}'
        logfc_suffix = '' if logfc_filter is None else f'_logFC{logfc_filter}'
        output_file = os.path.join(output_dir, f"{base_name}{position_suffix}{length_suffix}{logfc_suffix}_logo.pdf")
        
        # 创建标题
        title_parts = []
        title_parts.append(base_name)
        if position != 'all':
            title_parts.append(f"{position}端")
        if length_filter is not None:
            title_parts.append(f"长度={length_filter}nt")
        if logfc_filter is not None:
            title_parts.append(f"logFC{logfc_filter}")
        title = ' '.join(title_parts)
        
        # 创建序列logo图
        create_sequence_logo(freq_matrix, output_file, title)
        
    except Exception as e:
        print(f"处理文件 {file_path} 时出错: {str(e)}")

def main():
    args = parse_args()
    
    # 获取输入目录中所有以_length.csv结尾的文件
    input_pattern = os.path.join(args.input_dir, '*_length.csv')
    csv_files = glob.glob(input_pattern)
    
    if not csv_files:
        print(f"错误: 在 {args.input_dir} 中没有找到以_length.csv结尾的文件")
        sys.exit(1)
    
    # 确保输出目录存在
    os.makedirs(args.input_dir, exist_ok=True)
    
    # 处理每个CSV文件
    for file_path in csv_files:
        process_csv_file(
            file_path, 
            args.input_dir, 
            args.length_filter, 
            args.top_n, 
            args.position,
            args.logfc_filter
        )
    
    print(f"已完成 {len(csv_files)} 个文件的处理")

if __name__ == '__main__':
    main()
