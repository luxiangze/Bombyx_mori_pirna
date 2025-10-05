#!/usr/bin/env python3
"""
从CSV文件重新绘制piRNA分析图表，避免重新解析map文件
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import argparse
import glob
from collections import defaultdict

# 设置图表样式为浅色系，与原始脚本一致
plt.style.use('seaborn-v0_8-pastel')

def extract_sample_name(file_path):
    """从文件路径中提取样本名称"""
    file_name = os.path.basename(file_path)
    # 移除常见后缀
    for suffix in ['_length_distribution.csv', '_5prime_shift.csv', '_3prime_position.csv', '_unique_reads_length.csv']:
        if file_name.endswith(suffix):
            file_name = file_name[:-len(suffix)]
            break
    # 移除其他常见后缀
    for suffix in ['.fa.collapsed', '.map', '.no-dust']:
        if file_name.endswith(suffix):
            file_name = file_name[:-len(suffix)]
    return file_name

def load_csv_data(input_dir):
    """加载指定目录下的所有CSV文件"""
    data = {}
    
    # 查找所有长度分布CSV文件
    length_files = glob.glob(os.path.join(input_dir, '*_length_distribution.csv'))
    shift_files = glob.glob(os.path.join(input_dir, '*_5prime_shift.csv'))
    position_files = glob.glob(os.path.join(input_dir, '*_3prime_position.csv'))
    unique_length_files = glob.glob(os.path.join(input_dir, '*_unique_reads_length.csv'))
    
    # 为每个样本加载数据
    for length_file in length_files:
        sample_name = extract_sample_name(length_file)
        
        # 查找对应的其他CSV文件
        shift_file = None
        for f in shift_files:
            if extract_sample_name(f) == sample_name:
                shift_file = f
                break
        
        position_file = None
        for f in position_files:
            if extract_sample_name(f) == sample_name:
                position_file = f
                break
        
        unique_length_file = None
        for f in unique_length_files:
            if extract_sample_name(f) == sample_name:
                unique_length_file = f
                break
        
        # 加载数据
        df_length = pd.read_csv(length_file)
        df_shift = pd.read_csv(shift_file) if shift_file else None
        df_position = pd.read_csv(position_file) if position_file else None
        df_unique_length = pd.read_csv(unique_length_file) if unique_length_file else None
        
        data[sample_name] = (df_length, df_shift, df_position, df_unique_length)
    
    return data

def plot_combined_distributions(all_data, output_prefix):
    """绘制组合分布图"""
    # 设置图表样式为浅色系，与原始脚本一致
    plt.style.use('seaborn-v0_8-pastel')
    
    # 浅色系颜色和样式 - 扩展以支持更多样本
    colors = ['#8dd3c7', '#fb8072', '#bebada', '#80b1d3', '#fdb462', '#b3de69', '#fccde5', '#d9d9d9', '#bc80bd', '#ccebc5']
    markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h']
    linestyles = ['-', '--', '-.', ':', '-', '--', '-.', ':', '-', '--']
    
    # 计算样本数量和柱状图宽度
    num_samples = len(all_data)
    bar_width = 0.8 / num_samples if num_samples > 1 else 0.6
    
    # 简化样本名称，去除多余后缀
    simplified_sample_names = {}
    for sample_name in all_data.keys():
        # 去除".fa.collapsed"等后缀
        simple_name = sample_name
        for suffix in ['.map', '.fa.collapsed.map', '.fa']:
            if simple_name.endswith(suffix):
                simple_name = simple_name[:-len(suffix)]
                break
        simplified_sample_names[sample_name] = simple_name
    
    # 1. 绘制长度分布图
    plt.figure(figsize=(10, 6))
    for i, (sample_name, (df_length, _, _, _)) in enumerate(all_data.items()):
        total_count = df_length['Count'].sum()
        percentage = df_length['Count'] / total_count * 100 if total_count > 0 else df_length['Count']
        x = df_length['Length'] + i * bar_width - (num_samples * bar_width / 2) + (bar_width / 2)
        plt.bar(x, percentage, width=bar_width, color=colors[i % len(colors)], 
                label=simplified_sample_names[sample_name])
    
    # 确保x轴刻度为整数
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(1))
    
    plt.xlabel('Length (nt)')
    plt.ylabel('Percentage (%)')
    plt.title('Length Distribution')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_prefix}length_distribution.pdf", format='pdf')
    plt.close()
    
    # 2. 绘制5'端位移分布图
    plt.figure(figsize=(10, 6))
    for i, (sample_name, (_, df_shift, _, _)) in enumerate(all_data.items()):
        if df_shift is None:
            continue
        # 只保留 -2, -1, 0, 1, 2 五个值
        filtered_df = df_shift[df_shift['Shift'].isin([-2, -1, 0, 1, 2])].copy()
        plt.plot(filtered_df['Shift'], filtered_df['Percentage'], 
                 marker=markers[i % len(markers)], linestyle=linestyles[i % len(linestyles)], 
                 color=colors[i % len(colors)], linewidth=2, 
                 label=simplified_sample_names[sample_name])
    plt.xlabel("5' position shift")
    plt.ylabel('% of analyzed piRNA loci')
    plt.title("5' End Shift")
    # 设置 x 轴仅显示 -2, -1, 0, 1, 2
    plt.xticks([-2, -1, 0, 1, 2])
    # 确保 y 轴从 0 开始
    plt.ylim(bottom=0)
    # 确保x轴刻度为整数
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(1))
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_prefix}5prime_shift.pdf", format='pdf')
    plt.close()
    
    # 3. 绘制 3' 端相对位置分布图
    plt.figure(figsize=(10, 6))
    for i, (sample_name, (_, _, df_position, _)) in enumerate(all_data.items()):
        if df_position is None:
            continue
        plt.plot(df_position['Position'], df_position['Percentage'], 
                 marker=markers[i % len(markers)], linestyle=linestyles[i % len(linestyles)], 
                 color=colors[i % len(colors)], linewidth=2, 
                 label=simplified_sample_names[sample_name])
    
    plt.xlabel("3' relative position")
    plt.ylabel('% of analyzed piRNA loci')
    plt.title("3' End Position")
    # 确保x轴刻度为整数
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(1))
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_prefix}3prime_position.pdf", format='pdf')
    plt.close()
    
    # 4. 绘制唯一读数长度分布图
    plt.figure(figsize=(10, 6))
    for i, (sample_name, (_, _, _, df_unique_length)) in enumerate(all_data.items()):
        if df_unique_length is None:
            continue
            
        # 简化样本名称，去除多余后缀
        simplified_name = simplified_sample_names[sample_name]
        # 移除测序相关的后缀
        patterns_to_remove = ['_1.fq', '_2.fq', '.fq', '_1.fastq', '_2.fastq', '.fastq']
        for pattern in patterns_to_remove:
            if pattern in simplified_name:
                simplified_name = simplified_name.replace(pattern, '')
        
        # 如果样本名称包含"-KD"，只保留该部分
        if '-KD' in simplified_name:
            parts = simplified_name.split('-KD')
            if len(parts) > 1:
                simplified_name = parts[0] + '-KD' + parts[1].split('_')[0]
                
        total_unique = df_unique_length['UniqueCount'].sum()
        percentage = df_unique_length['UniqueCount'] / total_unique * 100 if total_unique > 0 else df_unique_length['UniqueCount']
        x = df_unique_length['Length'] + i * bar_width - (num_samples * bar_width / 2) + (bar_width / 2)
        plt.bar(x, percentage, width=bar_width, color=colors[i % len(colors)], 
                label=simplified_name)
    
    plt.xlabel('Length (nt)')
    plt.ylabel('Percentage (%)')
    plt.title('Unique Reads Length Distribution')
    # 确保x轴刻度为整数
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(1))
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_prefix}unique_length_distribution.pdf", format='pdf')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='从CSV文件重新绘制piRNA分析图表')
    parser.add_argument('-i', '--input', required=True, help='包含CSV文件的输入目录')
    parser.add_argument('-o', '--output', required=True, help='输出目录')
    parser.add_argument('-p', '--prefix', default='', help='输出文件前缀')
    args = parser.parse_args()
    
    # 确保输出目录存在
    os.makedirs(args.output, exist_ok=True)
    
    # 加载CSV数据
    all_data = load_csv_data(args.input)
    
    if not all_data:
        print("未找到有效的CSV文件！")
        return
    
    # 设置输出前缀
    output_prefix = os.path.join(args.output, args.prefix)
    
    # 绘制组合分布图
    plot_combined_distributions(all_data, output_prefix)
    
    print(f"已成功重新绘制图表并保存到 {args.output} 目录")

if __name__ == "__main__":
    main()
