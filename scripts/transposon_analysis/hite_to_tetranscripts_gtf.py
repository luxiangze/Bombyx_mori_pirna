#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
将HiTE预测的转座子转换成符合TEtranscript要求的gtf格式
用法: python hite_to_tetranscripts_gtf.py <输入gff文件> <输出gtf文件>
"""

import sys
import os

def parse_attributes(attr_str):
    """解析GFF属性字段"""
    attrs = {}
    for field in attr_str.split(';'):
        if not field.strip():
            continue
        key, value = field.split('=', 1)
        attrs[key] = value
    return attrs

def format_gtf_attributes(attrs):
    """格式化GTF属性字段"""
    result = []
    for key, value in attrs.items():
        result.append(f'{key} "{value}"')
    return '; '.join(result)

def convert_hite_to_tetranscripts(input_file, output_file):
    """将HiTE GFF文件转换为TEtranscript GTF格式"""
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            line = line.strip()
            # 跳过注释行
            if line.startswith('#'):
                continue

            # 解析GFF行
            fields = line.split('\t')
            if len(fields) != 9:
                continue

            chrom, source, feature, start, end, score, strand, phase, attrs_str = fields

            # 解析属性
            attrs = parse_attributes(attrs_str)

            # 获取必要的属性
            te_id = attrs.get('id', '')
            te_name = attrs.get('name', '')
            classification = attrs.get('classification', 'Unknown')

            # 处理分类信息
            if '/' in classification:
                class_id, family_id = classification.split('/', 1)
            else:
                # 处理Unknown情况
                class_id = classification
                family_id = "Unknown"

            # 构建新的GTF属性
            new_attrs = {
                "gene_id": te_name,
                "transcript_id": f"{te_id}",
                "class_id": class_id,
                "family_id": family_id
            }

            # 格式化GTF行
            new_line = [
                chrom,
                "HiTE",  # 修改来源标识
                "exon",       # 所有特征都改为exon
                start,
                end,
                score,
                strand,
                phase,
                format_gtf_attributes(new_attrs)
            ]

            outfile.write('\t'.join(new_line) + '\n')

    print(f"转换完成！输出文件: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"用法: {sys.argv[0]} <输入gff文件> <输出gtf文件>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    if not os.path.exists(input_file):
        print(f"错误: 输入文件 '{input_file}' 不存在!")
        sys.exit(1)

    convert_hite_to_tetranscripts(input_file, output_file)
