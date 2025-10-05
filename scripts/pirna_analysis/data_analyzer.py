#!/usr/bin/env python3
"""
数据分析模块 - 负责对解析后的数据进行各种分析
"""

import pandas as pd
from collections import defaultdict

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
    counts = [length_counts[l] for l in lengths]
    
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
        for pirna_id, start_pos in read_to_pirnas[read_key]:
            shift = start_pos
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
    """分析3'端相对位置分布"""
    # 存储每个3'端相对位置的归一化计数
    position_counts = defaultdict(float)
    
    # 对每条读取序列进行处理
    for read_key, (count, read_length) in read_info.items():
        # 获取该读取序列比对到的参考piRNA数量
        pirna_count = read_to_pirna_count[read_key]
        # 归一化权重：读取序列计数 / 比对到的参考piRNA数量
        normalized_weight = count / pirna_count
        
        # 更新每个3'端相对位置的归一化计数
        for pirna_id, start_pos in read_to_pirnas[read_key]:
            # 计算3'端位置：起始位置 + 读取长度
            # 在seqmap格式中，trans_coord已经是从1开始的位置
            end_pos = start_pos + read_length
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

def analyze_base_composition_by_position(read_info, read_to_pirna_count, read_to_pirnas, base_range=(-2, 46)):
    """
    统计每个相对位置的碱基分布百分比。
    - base_range: (start, end)，统计的相对位置范围（包含端点）
    返回：
      碱基分布DataFrame（行：相对位置，列：A/T/C/G/N等百分比）
    """
    from collections import Counter
    base_pos_dict = {pos: Counter() for pos in range(base_range[0], base_range[1]+1)}
    base_total = {pos: 0 for pos in range(base_range[0], base_range[1]+1)}
    base_types = set()

    # 统计每个位点的碱基分布
    for read_key, (count, read_length) in read_info.items():
        pirna_count = read_to_pirna_count[read_key]
        normalized_weight = count / pirna_count if pirna_count > 0 else 0
        seq = read_key
        for rel_pos in range(base_range[0], base_range[1]+1):
            seq_pos = rel_pos if rel_pos >= 0 else rel_pos + len(seq)
            if 0 <= seq_pos < len(seq):
                base = seq[seq_pos].upper()
                base_pos_dict[rel_pos][base] += normalized_weight
                base_total[rel_pos] += normalized_weight
                base_types.add(base)

    # 生成碱基分布百分比DataFrame
    base_types = sorted(list(base_types))
    base_dist_rows = []
    for rel_pos in range(base_range[0], base_range[1]+1):
        total = base_total[rel_pos]
        row = {'Position': rel_pos}
        for base in base_types:
            pct = (base_pos_dict[rel_pos][base] / total * 100) if total > 0 else 0
            row[base] = pct
        base_dist_rows.append(row)
    df_base_dist = pd.DataFrame(base_dist_rows)

    return df_base_dist

    
    
    
