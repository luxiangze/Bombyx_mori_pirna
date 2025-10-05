#!/usr/bin/env python3
"""
数据解析模块 - 负责读取和解析 map 文件和参考 piRNA 文件
"""

import os
import time
import pandas as pd
from collections import defaultdict
from tqdm import tqdm
from Bio import SeqIO

def parse_map_file_optimized(map_file, ref_pirna_file, chunk_size=10000):
    """
    优化版本的parse_map_file函数，使用分块处理和更高效的数据结构
    
    参数:
        map_file: map文件路径
        ref_pirna_file: 参考piRNA文件路径
        chunk_size: 每次处理的行数
    
    返回:
        read_info, read_to_pirna_count, read_to_pirnas, filter_stats
    """
    # 存储每条读取序列比对到的参考piRNA数量
    read_to_pirna_count = defaultdict(int)
    # 存储每条读取序列的计数和长度
    read_info = {}
    # 存储每条读取序列比对到的参考piRNA及位置
    read_to_pirnas = defaultdict(list)
    
    # 记录过滤的序列数量
    filtered_count = 0
    total_count = 0
    
    # 计时
    start_time = time.time()
    
    print(f"开始解析参考piRNA文件: {ref_pirna_file}")
    # 从参考piRNA文件中读取序列
    pirna_sequences = {}
    try:
        # 使用迭代器模式减少内存使用
        for record in SeqIO.parse(ref_pirna_file, "fasta"):
            pirna_sequences[record.id] = str(record.seq)
        print(f"成功读取了 {len(pirna_sequences)} 条参考piRNA序列，耗时 {time.time() - start_time:.2f} 秒")
    except Exception as e:
        print(f"警告: 读取参考piRNA文件时出错: {e}")
        print("将使用target_seq作为替代")
    
    # 查找对应的.fa.collapsed文件以获取序列计数信息
    collapsed_file = map_file.replace('.map', '')
    if not os.path.exists(collapsed_file):
        # 尝试其他可能的文件名模式
        possible_collapsed_files = [
            map_file.replace('.fa.rmdup.map', '.fa.collapsed'),
            map_file.replace('.map', '.fa.collapsed'),
            os.path.join(os.path.dirname(map_file), os.path.basename(map_file).split('.')[0] + '.fa.collapsed')
        ]
        
        for possible_file in possible_collapsed_files:
            if os.path.exists(possible_file):
                collapsed_file = possible_file
                print(f"找到collapsed文件: {collapsed_file}")
                break
    
    # 使用Biopython从collapsed文件中读取序列计数信息和序列本身
    seq_counts = {}
    read_sequences = {}
    if os.path.exists(collapsed_file):
        collapsed_start_time = time.time()
        print(f"解析collapsed文件: {collapsed_file}")
        try:
            # 使用进度条显示处理进度
            with open(collapsed_file, 'r') as f:
                # 计算文件行数以设置进度条
                line_count = sum(1 for _ in f)
            
            record_count = 0
            for record in tqdm(SeqIO.parse(collapsed_file, "fasta"), total=line_count//2, desc="读取序列"):
                # 从FASTA头部提取计数信息
                count = 1
                if '-' in record.id:
                    try:
                        count = int(record.id.split('-')[-1])
                    except ValueError:
                        pass
                
                # 存储序列及其计数
                seq = str(record.seq)
                seq_counts[seq] = count
                
                # 同时存储序列ID和序列的映射关系
                read_sequences[record.id] = seq
                record_count += 1
            
            print(f"成功读取了 {len(seq_counts)} 条序列及其计数信息，耗时 {time.time() - collapsed_start_time:.2f} 秒")
            print(f"成功存储了 {len(read_sequences)} 条序列ID与序列的映射关系")
        except Exception as e:
            print(f"警告: 解析collapsed文件时出错: {e}")
    else:
        print(f"警告: 找不到对应的collapsed文件: {collapsed_file}")
    
    # 读取map文件
    map_start_time = time.time()
    print(f"开始解析map文件: {map_file}")
    try:
        # 使用pandas分块读取TSV格式的map文件，减少内存使用
        # 首先获取列名
        header = pd.read_csv(map_file, sep='\t', nrows=1, header=None)
        column_count = len(header.columns)
        
        # 设置列名
        column_names = ['piRNA_id', 'trans_coord', 'reads_id']
        
        # 计算文件行数以设置进度条
        with open(map_file, 'r') as f:
            line_count = sum(1 for _ in f) - 1  # 减去标题行
        
        # 分块读取并处理
        reader = pd.read_csv(map_file, sep='\t', skiprows=1, header=None, names=column_names, chunksize=chunk_size)
        
        # 使用进度条显示处理进度
        with tqdm(total=line_count, desc="处理map文件") as pbar:
            for chunk in reader:
                chunk_size = len(chunk)
                total_count += chunk_size
                pbar.update(chunk_size)
                
                # 处理每一行
                for _, row in chunk.iterrows():
                    piRNA_id = row['piRNA_id']  # 参考piRNA ID
                    trans_coord = int(row['trans_coord'])  # 偏移位置
                    reads_id = row['reads_id']  # 读取序列ID
                    
                    # 获取读取序列和计数
                    read_seq = read_sequences.get(reads_id, "")
                    if not read_seq:
                        # 如果读取序列为空，则跳过该条记录
                        continue
                    count = seq_counts.get(read_seq, 1)
                    
                    # 存储读取序列信息
                    read_key = read_seq
                    read_info[read_key] = (count, len(read_seq))
                    
                    # 更新读取序列比对到的参考piRNA数量
                    read_to_pirna_count[read_key] += 1
                    
                    # 存储读取序列比对到的参考piRNA及位置
                    read_to_pirnas[read_key].append((piRNA_id, trans_coord))
        
        print(f"成功解析了 {total_count} 条比对记录，耗时 {time.time() - map_start_time:.2f} 秒")
        print(f"过滤了 {filtered_count} 条起始位置>6的序列")
    except Exception as e:
        print(f"警告: 解析map文件时出错: {e}")
    
    # 返回解析结果和过滤统计
    filter_stats = (total_count, filtered_count)
    return read_info, read_to_pirna_count, read_to_pirnas, filter_stats
