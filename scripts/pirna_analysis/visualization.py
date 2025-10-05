#!/usr/bin/env python3
"""
可视化模块 - 负责生成各种图表
"""

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import os
import pandas as pd
import logomaker

def plot_combined_distributions(all_data, output_prefix):
    """绘制多个样本的分布图，将三个子图分别输出为单独的图片"""
    # 设置图表样式为浅色系
    plt.style.use('seaborn-v0_8-pastel')
    
    # 浅色系颜色和样式 - 扩展以支持更多样本
    colors = ['#8dd3c7', '#fb8072', '#bebada', '#80b1d3', '#fdb462', '#b3de69', '#fccde5', '#d9d9d9', '#bc80bd', '#ccebc5']
    linestyles = ['-', '--', '-.', ':', '-', '--', '-.', ':', '-', '--']
    markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h']
    
    # 计算最佳柱形图宽度
    num_samples = len(all_data)
    bar_width = 0.8 / num_samples if num_samples > 0 else 0.25
    
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
    
    # 确保输出目录存在
    os.makedirs(os.path.dirname(output_prefix) if os.path.dirname(output_prefix) else output_prefix, exist_ok=True)
    
    # 1. 绘制长度分布图（纵坐标为百分比）
    plt.figure(figsize=(10, 6))
    for i, (sample_name, (df_length, _, _)) in enumerate(all_data.items()):
        total_count = df_length['Count'].sum()
        percent = (df_length['Count'] / total_count) * 100 if total_count > 0 else df_length['Count']
        x = df_length['Length'] + i * bar_width - (num_samples * bar_width / 2) + (bar_width / 2)
        plt.bar(x, percent, width=bar_width, color=colors[i % len(colors)], 
                label=simplified_sample_names[sample_name])
    
    # 确保x轴刻度为整数
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(1))
    
    plt.xlabel('Length (nt)')
    plt.ylabel('Percentage of reads (%)')
    plt.title('Length Distribution')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_prefix}length_distribution.pdf", format='pdf')
    plt.close()
    
    # 2. 绘制 5' 端位移分布图
    plt.figure(figsize=(10, 6))
    for i, (sample_name, (_, df_shift, _)) in enumerate(all_data.items()):
        # 只保留 0，1，2 三个值
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
    for i, (sample_name, (_, _, df_position)) in enumerate(all_data.items()):
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

def plot_base_composition_logo(df_base, output_file, title=None):
    """
    使用logomaker绘制碱基组成序列logo图。
    df_base: 行为位置，列为A/T/C/G/N等百分比（0-100），数值可为百分比或小数
    output_file: 输出图片路径（建议pdf或png）
    title: 可选标题
    """
    # 保证Position为index，只保留A/C/G/U/T列
    if 'Position' in df_base.columns:
        df_plot = df_base.set_index('Position')
    else:
        df_plot = df_base.copy()
    base_cols = [col for col in df_plot.columns if col in ('A','C','G','U','T')]
    df_plot = df_plot[base_cols]
    # 若为百分比，转为小数
    if df_plot.values.max() > 1.1:
        df_plot = df_plot / 100
    plt.figure(figsize=(min(16, 0.4*len(df_plot)), 4))
    logomaker.Logo(df_plot, shade_below=.5, fade_below=.5, font_name='Arial')
    plt.xlabel('Position')
    plt.ylabel('Frequency')
    if title:
        plt.title(title)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def plot_combined_unique_reads_length(all_unique_data, output_prefix):
    """绘制多个样本的去重后reads长度分布组合图"""
    # 设置图表样式为浅色系
    plt.style.use('seaborn-v0_8-pastel')
    
    # 创建图表
    plt.figure(figsize=(12, 7))
    
    # 浅色系颜色和样式 - 扩展以支持更多样本
    colors = ['#8dd3c7', '#fb8072', '#bebada', '#80b1d3', '#fdb462', '#b3de69', '#fccde5', '#d9d9d9', '#bc80bd', '#ccebc5']
    
    # 计算最佳柱形图宽度
    num_samples = len(all_unique_data)
    bar_width = 0.8 / num_samples if num_samples > 0 else 0.25
    
    # 绘制长度分布（纵坐标为百分比）
    for i, (sample_name, df_unique_length) in enumerate(all_unique_data.items()):
        # 简化样本名称，去除多余后缀
        simplified_name = sample_name
        # 移除常见的后缀
        for suffix in ['.map', '.fa.rmdup.map', '.fa.collapsed.map', '.fa.collapsed.no-dust.map', '.collapsed.no-dust.map', '.no-dust.map', '.fa.collapsed.no-dust', '.collapsed.no-dust', '.no-dust', '.fa']:
            if simplified_name.endswith(suffix):
                simplified_name = simplified_name[:-len(suffix)]
                break
        
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
        percent_unique = (df_unique_length['UniqueCount'] / total_unique) * 100 if total_unique > 0 else df_unique_length['UniqueCount']
        x = df_unique_length['Length'] + i * bar_width - (num_samples * bar_width / 2) + (bar_width / 2)
        plt.bar(x, percent_unique, width=bar_width, color=colors[i % len(colors)], label=simplified_name)
    
    plt.xlabel('Length (nt)')
    plt.ylabel('Percentage of unique sequences (%)')
    plt.title('Unique Reads Length Distribution Comparison')
    # 确保x轴刻度为整数
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(1))
    plt.legend()
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图表，使用简洁的文件名，PDF格式
    plt.savefig(f"{output_prefix}unique_length_distribution.pdf", format='pdf')
    plt.close()
