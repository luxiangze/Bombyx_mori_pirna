#!/usr/bin/env python3
"""
工具函数模块 - 提供各种通用的工具函数
"""

import os

def extract_sample_name(file_path):
    """从文件路径中提取简洁的样本名称"""
    # 获取文件名（不含路径）
    file_name = os.path.basename(file_path)
    
    # 移除文件扩展名和常见的后缀
    sample_name = file_name
    for suffix in ['.fa.rmdup.map', '.fa.collapsed.map', '.fa.collapsed.no-dust.map', '.collapsed.no-dust.map', '.no-dust.map', '.fa.collapsed.no-dust', '.collapsed.no-dust', '.no-dust', '.fa''.map']:
        if sample_name.endswith(suffix):
            sample_name = sample_name[:-len(suffix)]
            break
    
    # 移除"_1.fq"和类似的后缀
    patterns_to_remove = ['_1.fq', '_2.fq', '.fq', '_1.fastq', '_2.fastq', '.fastq']
    for pattern in patterns_to_remove:
        if pattern in sample_name:
            sample_name = sample_name.replace(pattern, '')
    
    # 如果样本名称包含"-KD"，只保留该部分
    if '-KD' in sample_name:
        parts = sample_name.split('-KD')
        if len(parts) > 1:
            sample_name = parts[0] + '-KD' + parts[1].split('_')[0]
            
    return sample_name
