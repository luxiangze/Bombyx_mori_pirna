#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 导入plot_pirna_trim函数的同时，创建一个修复版本
from pirna_analysis.visualization import plot_pirna_trim as original_plot_pirna_trim

# 创建一个修复版本的plot_pirna_trim函数
def fixed_plot_pirna_trim(all_pirna_trim_data, output_prefix):
    """修复版本的plot_pirna_trim函数，修复了颜色条相关的错误"""
    # 设置图表样式为浅色系
    plt.style.use('seaborn-v0_8-pastel')
    
    # 获取样本列表
    sample_names = list(all_pirna_trim_data.keys())
    n_samples = len(sample_names)
    
    # 如果样本数量小于2，无法进行比较
    if n_samples < 2:
        print("警告：至少需要2个样本才能进行piRNA剪切情况的比较分析")
        return
    
    # 简化样本名称，去除多余后缀
    simplified_sample_names = {}
    for sample_name in sample_names:
        # 去除常见的后缀
        simple_name = sample_name
        for suffix in ['.map', '.fa.rmdup.map', '.fa.collapsed.map', '.fa.collapsed.no-dust.map', 
                      '.collapsed.no-dust.map', '.no-dust.map', '.fa.collapsed.no-dust', 
                      '.collapsed.no-dust', '.no-dust', '.fa']:
            if simple_name.endswith(suffix):
                simple_name = simple_name[:-len(suffix)]
                break
        
        # 移除测序相关的后缀
        patterns_to_remove = ['_1.fq', '_2.fq', '.fq', '_1.fastq', '_2.fastq', '.fastq']
        for pattern in patterns_to_remove:
            if pattern in simple_name:
                simple_name = simple_name.replace(pattern, '')
        
        # 如果样本名称包含"-KD"，只保留该部分
        if '-KD' in simple_name:
            parts = simple_name.split('-KD')
            if len(parts) > 1:
                simple_name = parts[0] + '-KD' + parts[1].split('_')[0]
        
        simplified_sample_names[sample_name] = simple_name
    
    # 颜色映射，用于表示fold change
    cmap = plt.cm.coolwarm
    
    # 绘制每对样本的比较，每对样本保存为单独的文件
    for i in range(n_samples):
        for j in range(i+1, n_samples):
            # 获取两个样本的简化名称
            sample1_name = simplified_sample_names[sample_names[i]]
            sample2_name = simplified_sample_names[sample_names[j]]
            
            # 创建单独的图表和轴对象
            fig, ax = plt.subplots(figsize=(6, 6), dpi=100)
            
            # 获取两个样本的数据
            df1 = all_pirna_trim_data[sample_names[i]]
            df2 = all_pirna_trim_data[sample_names[j]]
            
            # 合并两个样本的数据，以piRNA_ID为键
            merged_df = pd.merge(df1, df2, on='piRNA_ID', how='inner', suffixes=('_1', '_2'))
            
            # 计算每个点的颜色（基于fold change）
            # 添加一个小的值避免除以零
            epsilon = 1e-10
            fold_changes = merged_df['Trim_index_2'] / (merged_df['Trim_index_1'] + epsilon)
            
            # 计算fold change的最小值和最大值，用于动态调整颜色映射范围
            # 对于对数变化，使用对称的范围更合适
            # 如果最小值小于1，则最大值应该是1/最小值
            # 如果最大值大于1，则最小值应该是1/最大值
            fc_min = fold_changes.min()
            fc_max = fold_changes.max()
            
            # 确保范围是对称的（相对于1）
            if fc_min < 1 and 1/fc_min > fc_max:
                vmax = 1/fc_min
                vmin = fc_min
            else:
                vmax = fc_max
                vmin = 1/fc_max
            
            # 限制范围，避免过大或过小
            vmax = min(vmax, 10)  # 最大不超过10倍
            vmin = max(vmin, 0.1)  # 最小不低于0.1倍
            
            # 创建规范化对象
            norm = plt.Normalize(vmin=vmin, vmax=vmax)
            
            # 计算每个点的颜色
            colors = [cmap(norm(fc)) for fc in fold_changes]
            
            # 绘制散点图
            sc = ax.scatter(merged_df['Trim_index_1'], merged_df['Trim_index_2'], 
                      c=colors, alpha=0.7, s=15, edgecolors='none')
            
            # 添加对角线（y=x线）
            max_val = max(merged_df['Trim_index_1'].max(), merged_df['Trim_index_2'].max())
            min_val = min(merged_df['Trim_index_1'].min(), merged_df['Trim_index_2'].min())
            # 确保最小值为1，使用对数刻度
            min_val = max(1, min_val)
            
            # 绘制对角线
            ax.plot([min_val, max_val], [min_val, max_val], 'k-', linewidth=1)
            
            # 设置坐标轴为对数刻度
            ax.set_xscale('log')
            ax.set_yscale('log')
            
            # 设置坐标轴范围
            ax.set_xlim(min_val, max_val*1.1)
            ax.set_ylim(min_val, max_val*1.1)
            
            # 设置坐标轴标签
            ax.set_xlabel(f'{sample1_name}\n(Normalized reads)')
            ax.set_ylabel(f'{sample2_name}\n(Normalized reads)')
            
            # 创建ScalarMappable对象用于颜色条
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            
            # 添加颜色条，使用fig.colorbar而不是plt.colorbar
            cbar = fig.colorbar(sm, ax=ax)
            cbar.set_label(f'Fold change ({vmin:.2f}-{vmax:.2f})')
            
            # 调整布局
            plt.tight_layout()
            
            # 保存单独的图表文件，使用样本名称作为文件名
            output_filename = f"{output_prefix}{sample1_name}_vs_{sample2_name}"
            plt.savefig(f"{output_filename}.pdf", format='pdf', bbox_inches='tight')
            plt.close()

def main():
    """测试plot_pirna_trim函数"""
    # 设置输入和输出目录
    input_dir = "/home/gyk/project/ld_pirna/results/pirna_analysis_20250405_164942"
    output_dir = "/home/gyk/project/ld_pirna/results/test_plot_pirna_trim"
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 读取所有样本的pirna_trim数据
    all_pirna_trim_data = {}
    sample_files = [f for f in os.listdir(input_dir) if f.endswith("_pirna_trim.csv")]
    
    print(f"找到以下样本文件:")
    for sample_file in sample_files:
        sample_name = sample_file.replace("_pirna_trim.csv", "")
        print(f"  - {sample_file}")
        
        # 读取CSV文件
        file_path = os.path.join(input_dir, sample_file)
        df = pd.read_csv(file_path)
        
        # 打印数据的基本信息
        print(f"    行数: {len(df)}")
        print(f"    列名: {', '.join(df.columns)}")
        
        # 存储数据
        all_pirna_trim_data[sample_name] = df
    
    print(f"\n总共找到 {len(all_pirna_trim_data)} 个样本")
    
    # 调用修复版本的plot_pirna_trim函数
    print(f"\n开始绘制piRNA剪切情况的比较图...")
    fixed_plot_pirna_trim(all_pirna_trim_data, f"{output_dir}/")
    
    print(f"\n绘图完成，结果保存在: {output_dir}")

if __name__ == "__main__":
    main()
