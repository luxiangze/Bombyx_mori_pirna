#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
将CSV文件转换为FASTA格式

根据CSV文件中的第1-3列（序列、control和treatment的reads数），
生成两个FASTA文件，每个记录的ID格式为"序号-reads数"，
文件名根据CSV文件名中的两个样本名称命名。
'''

import os
import sys
import csv
import re
import argparse

def parse_args():
    '''
    解析命令行参数
    '''
    parser = argparse.ArgumentParser(description='将CSV文件转换为FASTA格式')
    parser.add_argument('csv_file', help='输入的CSV文件路径')
    parser.add_argument('-o', '--output_dir', 
                      default='/home/gyk/project/ld_pirna/data/bm_smRNA_pirna_diffrencial', 
                      help='输出目录，默认为/home/gyk/project/ld_pirna/data/bm_smRNA_pirna_diffrencial')
    return parser.parse_args()

def extract_sample_names(filename):
    '''
    从CSV文件名中提取样本名称
    例如：从 "Control_vs_SUGP1-KD2_edger_results_with_counts_FDR_0.05_length.csv" 中提取 "Control" 和 "SUGP1-KD2"
    '''
    base_name = os.path.basename(filename)
    match = re.search(r'(.+?)_vs_(.+?)_', base_name)
    if match:
        return match.group(1), match.group(2)
    else:
        # 如果无法从文件名中提取样本名称，使用默认名称
        return 'sample1', 'sample2'

def csv_to_fasta(csv_file, output_dir='/home/gyk/project/ld_pirna/data/bm_smRNA_pirna_diffrencial'):
    '''
    将CSV文件转换为FASTA格式
    '''
    # 从文件名中提取样本名称
    sample1, sample2 = extract_sample_names(csv_file)
    
    # 从文件名中提取样本名称的部分，例如 "Control_vs_SUGP1-KD2"
    base_name = os.path.basename(csv_file)
    match = re.search(r'(.+?)_edger', base_name)
    if match:
        name_prefix = match.group(1)
    else:
        # 如果无法提取，使用样本名称
        name_prefix = f"{sample1}_vs_{sample2}"
    
    # 创建输出目录（如果不存在）
    os.makedirs(output_dir, exist_ok=True)
    
    # 定义输出文件路径，使用简洁的文件名
    fasta_file1 = os.path.join(output_dir, f"{name_prefix}-control.fa.collapsed")
    fasta_file2 = os.path.join(output_dir, f"{name_prefix}-treatment.fa.collapsed")
    
    # 打开输出文件
    with open(fasta_file1, 'w') as f1, open(fasta_file2, 'w') as f2:
        # 读取CSV文件
        with open(csv_file, 'r') as csvfile:
            reader = csv.reader(csvfile)
            # 跳过标题行
            header = next(reader)
            
            # 检查CSV格式是否符合预期
            if len(header) < 3:
                print(f"错误：CSV文件格式不正确，至少需要3列数据，但只有{len(header)}列", file=sys.stderr)
                sys.exit(1)
            
            # 处理每一行数据
            for i, row in enumerate(reader, 1):
                if len(row) < 3:
                    print(f"警告：第{i}行数据不完整，跳过", file=sys.stderr)
                    continue
                
                # 提取序列和reads数
                sequence = row[0]
                control_reads = int(float(row[1])) if row[1] else 0
                treatment_reads = int(float(row[2])) if row[2] else 0
                
                # 写入FASTA格式
                if control_reads > 0:
                    f1.write(f">{i}-{control_reads}\n{sequence}\n")
                if treatment_reads > 0:
                    f2.write(f">{i}-{treatment_reads}\n{sequence}\n")
    
    print(f"已生成FASTA文件：\n{fasta_file1}\n{fasta_file2}")
    return fasta_file1, fasta_file2

def main():
    args = parse_args()
    csv_to_fasta(args.csv_file, args.output_dir)

if __name__ == '__main__':
    main()
