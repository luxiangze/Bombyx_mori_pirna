#!/usr/bin/env python3
"""
Analyze piRNA length distribution, 5' end shift distribution, and 3' end relative position.
Usage:
    python analyze_pirna_distribution.py -i <map_file1> [<map_file2> ...] -o <output_directory>
or:
    python analyze_pirna_distribution.py --input <map_file1> [<map_file2> ...] --output <output_directory>
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import time
import multiprocessing as mp
from collections import defaultdict
from Bio import SeqIO
from tqdm import tqdm

def parse_map_file(map_file, ref_pirna_file='/home/gyk/project/ld_pirna/data/bmo.v3.0.fa'):
    """Parse seqmap-generated map file and extract alignments to reference piRNAs.
    Filter out reads with alignment start position > 6 (original trans_coord > 9).

    Map file example (new format):
    trans_id        trans_coord     target_seq      probe_id        probe_seq       num_mismatch
    piR-bmo-1       1               TCCCGTATGGTCTAGTGGCTA   208883-4    TCCCGTATGGTCTAGTGGCTA   0

    Notes:
    1) We align positions 3–23 of each read; therefore start_pos = trans_coord - 3.
    2) ref_seq should be obtained from ref_pirna_file, not from target_seq.
    3) read_seq is in the collapsed file; use probe_id to retrieve it when available.
    """
    # 存储每条读取序列比对到的参考piRNA数量
    read_to_pirna_count = defaultdict(int)
    # 存储每条读取序列的计数和长度
    read_info = {}  # {read_key: (count, length)}
    # 存储每条读取序列比对到的参考piRNA及位置
    read_to_pirnas = defaultdict(list)
    # 存储读取序列ID和序列的映射关系
    read_sequences = {}
    
    # 记录过滤的序列数量
    filtered_count = 0
    total_count = 0
    
    print(f"Parsing reference piRNA file: {ref_pirna_file}")
    # 从参考piRNA文件中读取序列
    pirna_sequences = {}
    try:
        for record in SeqIO.parse(ref_pirna_file, "fasta"):
            pirna_sequences[record.id] = str(record.seq)
        print(f"Loaded {len(pirna_sequences)} reference piRNA sequences")
    except Exception as e:
        print(f"Warning: error reading reference piRNA file: {e}")
        print("Fallback to target_seq")
    
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
                print(f"Found collapsed file: {collapsed_file}")
                break
    
    # 使用Biopython从collapsed文件中读取序列计数信息和序列本身
    seq_counts = {}
    if os.path.exists(collapsed_file):
        print(f"Parsing collapsed file: {collapsed_file}")
        try:
            for record in SeqIO.parse(collapsed_file, "fasta"):
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
            
            print(f"Loaded {len(seq_counts)} sequences with counts")
            print(f"Stored {len(read_sequences)} read ID-to-sequence mappings")
        except Exception as e:
            print(f"Warning: error parsing collapsed file: {e}")
    else:
        print(f"Warning: collapsed file not found: {collapsed_file}")
    
    # 读取map文件
    print(f"Parsing map file: {map_file}")
    try:
        # 使用pandas读取TSV格式的map文件
        df_map = pd.read_csv(map_file, sep='\t', skiprows=1, header=None)
        # 确保至少有6列
        if df_map.shape[1] < 6:
            raise ValueError(f"Map file format error: fewer than 6 columns ({df_map.shape[1]})")
        
        # 重命名列以便更容易理解
        df_map.columns = ['trans_id', 'trans_coord', 'target_seq', 'probe_id', 'probe_seq', 'num_mismatch'] + \
                          [f'col{i+7}' for i in range(df_map.shape[1]-6)]
        
        # 处理每一行
        total_count = len(df_map)
        for _, row in df_map.iterrows():
            pirna_id = row['trans_id']  # 使用trans_id作为pirna_id
            raw_start_pos = int(row['trans_coord'])
            
            # 调整start_pos，因为使用每条序列的3-23nt的序列比对
            start_pos = raw_start_pos - 3
            
            # 从参考piRNA文件中获取ref_seq
            if pirna_id in pirna_sequences:
                ref_seq = pirna_sequences[pirna_id]
            else:
                # 如果找不到，使用target_seq作为替代
                ref_seq = row['target_seq']
                print(f"Warning: reference piRNA ID not found: {pirna_id}; using target_seq as fallback")
            
            probe_id = row['probe_id']
            read_seq = row['probe_seq']
            
            # 如果probe_id在read_sequences中有映射，优先使用映射的序列
            if probe_id in read_sequences:
                read_seq = read_sequences[probe_id]
            
            _mismatches = int(row['num_mismatch'])
            
            # 过滤比对起始位置大于6的序列（原始位置大于9）
            if raw_start_pos > 9:  # 相当于调整后的start_pos > 6
                filtered_count += 1
                continue
            
            # 从probe_id中提取计数信息，格式为ID-COUNT
            read_count = 1  # 默认计数为1
            if '-' in probe_id:
                try:
                    read_count = int(probe_id.split('-')[-1])
                except ValueError:
                    pass
            
            # 从seq_counts中获取序列计数（如果可用）
            if read_seq in seq_counts:
                read_count = seq_counts[read_seq]
            
            # 使用读取序列作为键，假定都是正向链
            strand = '+'
            read_key = f"{read_seq}_{strand}"
            
            # 更新读取序列的计数和长度
            read_info[read_key] = (read_count, len(read_seq))
            
            # 更新读取序列比对到的参考piRNA
            read_to_pirna_count[read_key] += 1
            read_to_pirnas[read_key].append((pirna_id, start_pos, len(ref_seq)))
    except Exception as e:
        print(f"Warning: error parsing map file: {e}")
    
    # 返回过滤统计信息以便在详细模式下显示
    filter_stats = (total_count, filtered_count)
    print(f"Parsing complete: total {total_count} records, filtered {filtered_count}")
    return read_info, read_to_pirna_count, read_to_pirnas, filter_stats

def analyze_length_distribution(read_info, read_to_pirna_count):
    """分析piRNA长度分布"""
    # 存储每个长度的总计数
    length_counts = defaultdict(int)
    
    # 对每条读取序列进行处理
    for read_key, (count, length) in read_info.items():
        # 获取该读取序列比对到的参考piRNA数量
        pirna_count = read_to_pirna_count[read_key]
        # 归一化权重：读取序列计数 / 比对到的参考piRNA数量
        normalized_weight = count / pirna_count
        
        # 更新长度分布
        length_counts[length] += normalized_weight
    
    # 转换为DataFrame
    lengths = sorted(length_counts.keys())
    counts = [length_counts[length_val] for length_val in lengths]
    
    df_length = pd.DataFrame({
        'Length': lengths,
        'Count': counts
    })
    
    return df_length

def analyze_5prime_shift(read_info, read_to_pirna_count, read_to_pirnas):
    """分析5'端位移分布
    
    在seqmap格式中，trans_coord表示比对起始位置（从1开始）
    """
    # 存储每个位移值的归一化计数
    shift_counts = defaultdict(float)
    
    # 对每条读取序列进行处理
    for read_key, (count, read_length) in read_info.items():
        # 获取该读取序列比对到的参考piRNA数量
        pirna_count = read_to_pirna_count[read_key]
        # 归一化权重：读取序列计数 / 比对到的参考piRNA数量
        normalized_weight = count / pirna_count
        
        # 更新每个位移值的归一化计数
        for pirna_id, start_pos, ref_length in read_to_pirnas[read_key]:
            # 5'端位移：起始位置 - 1 (因为位置从1开始)
            # 在seqmap格式中，trans_coord已经是从1开始的位置
            shift = start_pos - 1
            shift_counts[shift] += normalized_weight
    
    # 转换为DataFrame
    shifts = sorted(shift_counts.keys())
    counts = [shift_counts[s] for s in shifts]
    
    # 计算总计数
    total_count = sum(counts)
    
    # 计算百分比
    percentages = [count / total_count * 100 for count in counts] if total_count > 0 else [0] * len(counts)
    
    df_shift = pd.DataFrame({
        'Shift': shifts,
        'Count': counts,
        'Percentage': percentages
    })
    
    return df_shift

def analyze_3prime_position(read_info, read_to_pirna_count, read_to_pirnas):
    """分析3'端相对位置分布
    
    在seqmap格式中，需要计算3'端相对位置：起始位置(trans_coord) + 读取长度 - 1
    """
    # 存储每个3'端相对位置的归一化计数
    position_counts = defaultdict(float)
    
    # 对每条读取序列进行处理
    for read_key, (count, read_length) in read_info.items():
        # 获取该读取序列比对到的参考piRNA数量
        pirna_count = read_to_pirna_count[read_key]
        # 归一化权重：读取序列计数 / 比对到的参考piRNA数量
        normalized_weight = count / pirna_count
        
        # 更新每个3'端相对位置的归一化计数
        for pirna_id, start_pos, ref_length in read_to_pirnas[read_key]:
            # 计算3'端位置：起始位置 + 读取长度 - 1
            # 在seqmap格式中，trans_coord已经是从1开始的位置
            end_pos = start_pos + read_length - 1
            position_counts[end_pos] += normalized_weight
    
    # 转换为DataFrame
    positions = sorted(position_counts.keys())
    counts = [position_counts[p] for p in positions]
    
    # 计算总计数
    total_count = sum(counts)
    
    # 计算百分比
    percentages = [count / total_count * 100 for count in counts] if total_count > 0 else [0] * len(counts)
    
    df_position = pd.DataFrame({
        'Position': positions,
        'Count': counts,
        'Percentage': percentages
    })
    
    return df_position

def plot_distributions(df_length, df_shift, df_position, sample_name, output_prefix):
    """Plot distribution figures"""
    # 设置图表样式为浅色系
    plt.style.use('seaborn-v0_8-pastel')
    
    # 创建一个包含三个子图的图表
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
    
    # 绘制长度分布 (类似于图B)
    ax1.bar(df_length['Length'], df_length['Count'] / 1e6, color='#1f77b4', width=0.8)
    ax1.set_xlabel('Length (nt)')
    ax1.set_ylabel('Reads per million (x 10⁶)')
    ax1.set_title(f'{sample_name} - Length Distribution')
    
    # 绘制5'端位移分布 (类似于图C左侧)
    ax2.plot(df_shift['Shift'], df_shift['Percentage'], 'o-', color='#d62728', linewidth=2)
    ax2.set_xlabel("5' position shift")
    ax2.set_ylabel('% of analyzed piRNA loci')
    ax2.set_title(f"{sample_name} - 5' End Shift")
    
    # 绘制3'端相对位置分布 (类似于图C右侧)
    ax3.plot(df_position['Position'], df_position['Percentage'], 'o-', color='#ff7f0e', linewidth=2)
    ax3.set_xlabel("3' relative position")
    ax3.set_ylabel('% of analyzed piRNA loci')
    ax3.set_title(f"{sample_name} - 3' End Position")
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图表为PDF格式
    plt.savefig(f"{output_prefix}_{sample_name}_distributions.pdf", format='pdf')
    plt.close()

def plot_combined_distributions(all_data, output_prefix):
    """Plot combined distributions for multiple samples; export three figures"""
    # 设置图表样式为浅色系
    plt.style.use('seaborn-v0_8-pastel')
    
    # 浅色系颜色和样式 - 扩展以支持更多样本
    colors = ['#8dd3c7', '#fb8072', '#bebada', '#80b1d3', '#fdb462', '#b3de69', '#fccde5', '#d9d9d9', '#bc80bd', '#ccebc5']#ccebc5']
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
    
    # 1. 绘制长度分布图
    plt.figure(figsize=(10, 6))
    for i, (sample_name, (df_length, _, _)) in enumerate(all_data.items()):
        x = df_length['Length'] + i * bar_width - (num_samples * bar_width / 2) + (bar_width / 2)
        plt.bar(x, df_length['Count'] / 1e6, width=bar_width, color=colors[i % len(colors)], 
                label=simplified_sample_names[sample_name])
    
    plt.xlabel('Length (nt)')
    plt.ylabel('Reads per million (x 10⁶)')
    plt.title('Length Distribution')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_prefix}/length_distribution.pdf", format='pdf')
    plt.close()
    
    # 2. 绘制 5' 端位移分布图
    plt.figure(figsize=(10, 6))
    for i, (sample_name, (_, df_shift, _)) in enumerate(all_data.items()):
        # 只保留 0，1，2 三个值
        filtered_df = df_shift[df_shift['Shift'].isin([0, 1, 2])].copy()
        plt.plot(filtered_df['Shift'], filtered_df['Percentage'], 
                 marker=markers[i % len(markers)], linestyle=linestyles[i % len(linestyles)], 
                 color=colors[i % len(colors)], linewidth=2, 
                 label=simplified_sample_names[sample_name])
    
    plt.xlabel("5' position shift")
    plt.ylabel('% of analyzed piRNA loci')
    plt.title("5' End Shift")
    # 设置 x 轴仅显示 0，1，2
    plt.xticks([0, 1, 2])
    # 确保 y 轴从 0 开始
    plt.ylim(bottom=0)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_prefix}/5prime_shift.pdf", format='pdf')
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
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_prefix}/3prime_position.pdf", format='pdf')
    plt.close()

def analyze_unique_reads_length(read_info):
    """分析去重后的reads长度分布"""
    # 使用pandas的Series更高效地处理长度计数
    lengths = [length for _, (_, length) in read_info.items()]
    length_series = pd.Series(lengths)
    
    # 计算每个长度的唯一序列数量
    unique_length_counts = length_series.value_counts().sort_index()
    
    # 转换为DataFrame
    df_unique_length = pd.DataFrame({
        'Length': unique_length_counts.index,
        'UniqueCount': unique_length_counts.values
    })
    
    return df_unique_length

def plot_unique_reads_length(df_unique_length, sample_name, output_prefix):
    """绘制去重后的reads长度分布图"""
    # 设置图表样式为浅色系
    plt.style.use('seaborn-v0_8-pastel')
    
    # 创建图表
    plt.figure(figsize=(10, 6))
    
    # 绘制长度分布
    plt.bar(df_unique_length['Length'], df_unique_length['UniqueCount'], color='#b3de69', width=0.8)
    plt.xlabel('Length (nt)')
    plt.ylabel('Number of unique sequences')
    plt.title(f'{sample_name} - Unique Reads Length Distribution')
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图表为PDF格式
    plt.savefig(f"{output_prefix}_{sample_name}_unique_length_distribution.pdf", format='pdf')
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
    
    # 绘制长度分布
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
        
        x = df_unique_length['Length'] + i * bar_width - (num_samples * bar_width / 2) + (bar_width / 2)
        plt.bar(x, df_unique_length['UniqueCount'], width=bar_width, color=colors[i % len(colors)], label=simplified_name)
    
    plt.xlabel('Length (nt)')
    plt.ylabel('Number of unique sequences')
    plt.title('Unique Reads Length Distribution Comparison')
    plt.legend()
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图表，使用简洁的文件名，PDF格式
    plt.savefig(f"{output_prefix}/unique_length_distribution.pdf", format='pdf')
    plt.close()

def extract_sample_name(file_path):
    """从文件路径中提取简洁的样本名称"""
    # 获取文件名（不含路径）
    file_name = os.path.basename(file_path)
    
    # 移除文件扩展名和常见的后缀
    sample_name = file_name
    for suffix in ['.map', '.fa.rmdup.map', '.fa.collapsed.map', '.fa.collapsed.no-dust.map', '.collapsed.no-dust.map', '.no-dust.map', '.fa.collapsed.no-dust', '.collapsed.no-dust', '.no-dust', '.fa']:
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

def parse_map_file_optimized(map_file, ref_pirna_file='/home/gyk/project/ld_pirna/data/bmo.v3.0.fa', chunk_size=10000):
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
        column_names = ['trans_id', 'trans_coord', 'target_seq', 'probe_id', 'probe_seq', 'num_mismatch']
        if column_count > 6:
            column_names += [f'col{i+7}' for i in range(column_count-6)]
        
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
                
                for _, row in chunk.iterrows():
                    pirna_id = row['trans_id']  # 使用trans_id作为pirna_id
                    raw_start_pos = int(row['trans_coord'])
                    
                    # 调整start_pos，因为使用每条序列的3-23nt的序列比对
                    start_pos = raw_start_pos - 3
                    
                    # 从参考piRNA文件中获取ref_seq
                    if pirna_id in pirna_sequences:
                        ref_seq = pirna_sequences[pirna_id]
                    else:
                        # 如果找不到，使用target_seq作为替代
                        ref_seq = row['target_seq']
                        if total_count < 10:  # 只显示前几条警告，避免刷屏
                            print(f"警告: 找不到piRNA ID: {pirna_id}的参考序列，使用target_seq作为替代")
                    
                    probe_id = row['probe_id']
                    read_seq = row['probe_seq']
                    
                    # 如果probe_id在read_sequences中有映射，优先使用映射的序列
                    if probe_id in read_sequences:
                        read_seq = read_sequences[probe_id]
                    
                    _mismatches = int(row['num_mismatch'])
                    
                    # 过滤比对起始位置大于6的序列（原始位置大于9）
                    if raw_start_pos > 9:  # 相当于调整后的start_pos > 6
                        filtered_count += 1
                        continue
                    
                    # 从probe_id中提取计数信息，格式为ID-COUNT
                    read_count = 1  # 默认计数为1
                    if '-' in probe_id:
                        try:
                            read_count = int(probe_id.split('-')[-1])
                        except ValueError:
                            pass
                    
                    # 从seq_counts中获取序列计数（如果可用）
                    if read_seq in seq_counts:
                        read_count = seq_counts[read_seq]
                    
                    # 使用读取序列作为键，假定都是正向链
                    strand = '+'
                    read_key = f"{read_seq}_{strand}"
                    
                    # 更新读取序列的计数和长度
                    read_info[read_key] = (read_count, len(read_seq))
                    
                    # 更新读取序列比对到的参考piRNA
                    read_to_pirna_count[read_key] += 1
                    read_to_pirnas[read_key].append((pirna_id, start_pos, len(ref_seq)))
    except Exception as e:
        print(f"Warning: error parsing map file: {e}")
        import traceback
        traceback.print_exc()
    
    # 返回过滤统计信息以便在详细模式下显示
    filter_stats = (total_count, filtered_count)
    print(f"Parsing complete: total {total_count} records, filtered {filtered_count}, elapsed {time.time() - start_time:.2f} s")
    return read_info, read_to_pirna_count, read_to_pirnas, filter_stats

def process_file(map_file, ref_file, output_prefix, verbose=False):
    """
    处理单个map文件的函数，用于多进程并行处理
    只负责数据处理和保存，不进行绘图
    """
    sample_name = extract_sample_name(map_file)
    if verbose:
        print(f"Processing sample {sample_name}: {map_file}")
    
    # 使用优化版本的解析函数
    read_info, read_to_pirna_count, read_to_pirnas, filter_stats = parse_map_file_optimized(map_file, ref_file)
    
    if verbose:
        total_reads, filtered_reads = filter_stats
        print(f"  Found {len(read_info)} unique reads in {sample_name}")
        print(f"  Filtered {filtered_reads}/{total_reads} reads with start position > 9 (adjusted > 6)")
    
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
    
    # 返回样本名称和处理结果
    return sample_name

def main():
    # Configure CLI
    parser = argparse.ArgumentParser(description="Analyze piRNA length distribution, 5' end shift, 3' end relative position, and unique reads length distribution")
    parser.add_argument('-i', '--input', nargs='+', required=True, help='One or more input .map files')
    parser.add_argument('-o', '--output', required=True, help='Output directory or prefix')
    parser.add_argument('-v', '--verbose', action='store_true', help='Show verbose processing information')
    parser.add_argument('-p', '--prefix', default='', help='Output filename prefix (default: empty)')
    parser.add_argument('--no-biopython', action='store_true', help='Do not use Biopython; fall back to legacy parsing')
    parser.add_argument('-r', '--ref', default='/home/gyk/project/ld_pirna/data/bmo.v3.0.fa', help='Reference piRNA FASTA file path')
    parser.add_argument('--threads', type=int, default=0, help='Number of worker processes (0 uses all available CPU cores)')
    
    # 解析命令行参数
    args = parser.parse_args()
    
    # 获取输入文件列表和输出前缀
    map_files = args.input
    output_prefix = args.output
    verbose = args.verbose
    file_prefix = args.prefix
    use_biopython = not args.no_biopython
    
    # 设置多进程的线程数
    threads = args.threads
    if threads <= 0:
        threads = mp.cpu_count()
    
    if use_biopython:
        try:
            import Bio
            if verbose:
                print(f"Using Biopython {Bio.__version__} for sequence analysis")
        except ImportError:
            print("Warning: Biopython not found; using legacy parsing")
            use_biopython = False
    
    # 如果输出目录不存在，则创建
    output_dir = os.path.dirname(output_prefix) if os.path.dirname(output_prefix) else output_prefix
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    print(f"Using {threads} CPU core(s) for processing")
    
    # 处理多个样本
    all_data = {}
    all_unique_data = {}
    start_time = time.time()
    
    # 使用多进程并行处理文件
    if len(map_files) > 1 and threads > 1:
        print(f"Processing {len(map_files)} files in parallel...")
        # 创建进程池
        pool = mp.Pool(min(threads, len(map_files)))
        
        # 准备处理函数的参数
        process_args = [(map_file, args.ref, output_prefix, verbose) for map_file in map_files]
        
        # 并行处理文件
        results = []
        for result in tqdm(pool.starmap(process_file, process_args), total=len(map_files), desc="处理文件"):
            results.append(result)
        
        # 关闭进程池
        pool.close()
        pool.join()
        
        print("数据处理完成，开始生成图表...")
        
        # 收集所有样本的数据用于组合图
        for map_file in map_files:
            sample_name = extract_sample_name(map_file)
            
            # 读取已保存的数据文件
            df_length = pd.read_csv(f"{output_prefix}/{sample_name}_length_distribution.csv")
            df_shift = pd.read_csv(f"{output_prefix}/{sample_name}_5prime_shift.csv")
            df_position = pd.read_csv(f"{output_prefix}/{sample_name}_3prime_position.csv")
            df_unique_length = pd.read_csv(f"{output_prefix}/{sample_name}_unique_reads_length.csv")
            
            # 在主进程中绘制图表
            if verbose:
                print(f"  为{sample_name}生成图表...")
                
            # 绘制各种分布图
            try:
                # 使用matplotlib绘制图表
                plt.figure(figsize=(10, 6))
                plt.bar(df_length['Length'], df_length['Count'], color='#1f77b4')
                plt.xlabel('Length (nt)')
                plt.ylabel('Count')
                plt.title(f'{sample_name} - Length Distribution')
                plt.tight_layout()
                plt.savefig(f"{output_prefix}/{sample_name}_length_distribution.png")
                plt.close()
                
                # 5'端位移分布图
                plt.figure(figsize=(10, 6))
                plt.bar(df_shift['Shift'], df_shift['Percentage'], color='#ff7f0e')
                plt.xlabel("5' position shift")
                plt.ylabel('% of analyzed piRNA loci')
                plt.title(f"{sample_name} - 5' End Shift")
                plt.tight_layout()
                plt.savefig(f"{output_prefix}/{sample_name}_5prime_shift.png")
                plt.close()
                
                # 3'端相对位置分布图
                plt.figure(figsize=(10, 6))
                plt.plot(df_position['Position'], df_position['Percentage'], 'o-', color='#2ca02c', linewidth=2)
                plt.xlabel("3' relative position")
                plt.ylabel('% of analyzed piRNA loci')
                plt.title(f"{sample_name} - 3' End Position")
                plt.tight_layout()
                plt.savefig(f"{output_prefix}/{sample_name}_3prime_position.png")
                plt.close()
                
                # 去重后的reads长度分布图
                plt.figure(figsize=(10, 6))
                plt.bar(df_unique_length['Length'], df_unique_length['UniqueCount'], color='#d62728')
                plt.xlabel('Length (nt)')
                plt.ylabel('Number of unique sequences')
                plt.title(f'{sample_name} - Unique Reads Length Distribution')
                plt.tight_layout()
                plt.savefig(f"{output_prefix}/{sample_name}_unique_reads_length.png")
                plt.close()
            except Exception as e:
                print(f"警告: 为{sample_name}生成图表时出错: {e}")
            
            # 添加前缀并重命名文件（如果需要）
            if file_prefix:
                os.rename(f"{output_prefix}/{sample_name}_length_distribution.csv", 
                          f"{output_prefix}/{file_prefix}{sample_name}_length_distribution.tsv")
                os.rename(f"{output_prefix}/{sample_name}_5prime_shift.csv", 
                          f"{output_prefix}/{file_prefix}{sample_name}_5prime_shift.tsv")
                os.rename(f"{output_prefix}/{sample_name}_3prime_position.csv", 
                          f"{output_prefix}/{file_prefix}{sample_name}_3prime_position.tsv")
                os.rename(f"{output_prefix}/{sample_name}_unique_reads_length.csv", 
                          f"{output_prefix}/{file_prefix}{sample_name}_unique_length_distribution.tsv")
                
                # 重命名图片文件
                os.rename(f"{output_prefix}/{sample_name}_length_distribution.png", 
                          f"{output_prefix}/{file_prefix}{sample_name}_length_distribution.png")
                os.rename(f"{output_prefix}/{sample_name}_5prime_shift.png", 
                          f"{output_prefix}/{file_prefix}{sample_name}_5prime_shift.png")
                os.rename(f"{output_prefix}/{sample_name}_3prime_position.png", 
                          f"{output_prefix}/{file_prefix}{sample_name}_3prime_position.png")
                os.rename(f"{output_prefix}/{sample_name}_unique_reads_length.png", 
                          f"{output_prefix}/{file_prefix}{sample_name}_unique_reads_length.png")
            
            # 存储数据用于组合图
            all_data[sample_name] = (df_length, df_shift, df_position)
            all_unique_data[sample_name] = df_unique_length
    else:
        # 单进程处理文件
        for map_file in map_files:
            # 从文件名提取样本名称
            sample_name = extract_sample_name(map_file)
            
            if verbose:
                print(f"处理样本 {sample_name}: {map_file}")
            else:
                print(f"处理样本 {sample_name}...")
            
            # 使用优化版本的解析函数
            read_info, read_to_pirna_count, read_to_pirnas, filter_stats = parse_map_file_optimized(map_file, args.ref)
            
            if verbose:
                total_reads, filtered_reads = filter_stats
                print(f"  在{sample_name}中找到{len(read_info)}个唯一reads")
                print(f"  过滤了{filtered_reads}/{total_reads}个起始位置>9的reads（调整后>6）")
            
            # 分析长度分布
            df_length = analyze_length_distribution(read_info, read_to_pirna_count)
            
            # 分析5'端位移分布
            df_shift = analyze_5prime_shift(read_info, read_to_pirna_count, read_to_pirnas)
            
            # 分析3'端相对位置分布
            df_position = analyze_3prime_position(read_info, read_to_pirna_count, read_to_pirnas)
            
            # 分析去重后的reads长度分布
            df_unique_length = analyze_unique_reads_length(read_info)
            
            # 保存数据，使用简洁的文件名
            os.makedirs(f"{output_prefix}", exist_ok=True)
            df_length.to_csv(f"{output_prefix}/{file_prefix}{sample_name}_length_distribution.tsv", sep='\t', index=False)
            df_shift.to_csv(f"{output_prefix}/{file_prefix}{sample_name}_5prime_shift.tsv", sep='\t', index=False)
            df_position.to_csv(f"{output_prefix}/{file_prefix}{sample_name}_3prime_position.tsv", sep='\t', index=False)
            df_unique_length.to_csv(f"{output_prefix}/{file_prefix}{sample_name}_unique_length_distribution.tsv", sep='\t', index=False)
            
            if verbose:
                print(f"  为{sample_name}生成图表...")
            
            # 绘制图表
            try:
                # 使用matplotlib绘制图表
                plt.figure(figsize=(10, 6))
                plt.bar(df_length['Length'], df_length['Count'], color='#1f77b4')
                plt.xlabel('Length (nt)')
                plt.ylabel('Count')
                plt.title(f'{sample_name} - Length Distribution')
                plt.tight_layout()
                plt.savefig(f"{output_prefix}/{file_prefix}{sample_name}_length_distribution.png")
                plt.close()
                
                # 5'端位移分布图
                plt.figure(figsize=(10, 6))
                plt.bar(df_shift['Shift'], df_shift['Percentage'], color='#ff7f0e')
                plt.xlabel("5' position shift")
                plt.ylabel('% of analyzed piRNA loci')
                plt.title(f"{sample_name} - 5' End Shift")
                plt.tight_layout()
                plt.savefig(f"{output_prefix}/{file_prefix}{sample_name}_5prime_shift.png")
                plt.close()
                
                # 3'端相对位置分布图
                plt.figure(figsize=(10, 6))
                plt.plot(df_position['Position'], df_position['Percentage'], 'o-', color='#2ca02c', linewidth=2)
                plt.xlabel("3' relative position")
                plt.ylabel('% of analyzed piRNA loci')
                plt.title(f"{sample_name} - 3' End Position")
                plt.tight_layout()
                plt.savefig(f"{output_prefix}/{file_prefix}{sample_name}_3prime_position.png")
                plt.close()
                
                # 去重后的reads长度分布图
                plt.figure(figsize=(10, 6))
                plt.bar(df_unique_length['Length'], df_unique_length['UniqueCount'], color='#d62728')
                plt.xlabel('Length (nt)')
                plt.ylabel('Number of unique sequences')
                plt.title(f'{sample_name} - Unique Reads Length Distribution')
                plt.tight_layout()
                plt.savefig(f"{output_prefix}/{file_prefix}{sample_name}_unique_reads_length.png")
                plt.close()
            except Exception as e:
                print(f"警告: 为{sample_name}生成图表时出错: {e}")
            
            # 存储数据用于组合图
            all_data[sample_name] = (df_length, df_shift, df_position)
            all_unique_data[sample_name] = df_unique_length
    
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
