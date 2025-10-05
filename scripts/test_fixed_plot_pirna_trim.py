#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 添加脚本目录到Python路径
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# 导入修复后的plot_pirna_trim函数
from scripts.pirna_analysis.visualization import plot_pirna_trim

def main():
    """测试修复后的plot_pirna_trim函数"""
    # 设置输入和输出目录
    input_dir = "/home/gyk/project/ld_pirna/results/pirna_analysis_20250405_164942"
    output_dir = "/home/gyk/project/ld_pirna/results/test_fixed_plot_pirna_trim"
    
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
    
    # 调用修复后的plot_pirna_trim函数
    print(f"\n开始绘制piRNA剪切情况的比较图...")
    plot_pirna_trim(all_pirna_trim_data, f"{output_dir}/")
    
    print(f"\n绘图完成，结果保存在: {output_dir}")

if __name__ == "__main__":
    main()
