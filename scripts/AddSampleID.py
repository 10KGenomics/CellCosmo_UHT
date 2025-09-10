# -*- coding: utf-8 -*-
"""
AddSampleID.py

功能：
    从Cell_Summary.xls（纯文本、Tab分隔）文件中读取细胞级统计数据，
    根据CB列的barcode，将每个细胞分配到指定的Sample组。
    barcode库由外部文件mapping.txt给出（Tab分隔），第一列为barcode序列，
    第二列包含索引信息（如bc2-1表示索引1），范围和样本名称通过命令行参数灵活指定。
"""

import argparse
import sys
import pandas as pd
import re


def parse_args():
    """解析命令行参数并返回参数对象"""
    parser = argparse.ArgumentParser(
        description='Step9.AddSampleID: 根据CB列的barcode给每行添加SampleID（可变组数）',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '--inputfile', '-i', required=True,
        help='输入的Cell_Summary.xls（Tab分隔）文件路径')
    parser.add_argument(
        '--barcodelist', '-b', required=True,
        help='包含barcode的mapping.txt文件，Tab分隔，第一列为barcode，第二列含索引信息')
    parser.add_argument(
        '--splitCB', '-c', required=True,
        help="索引规范列表（1-based），支持：\n"
             "1. 范围格式：1-64（包含起始和结束）\n"
             "2. 单个索引：3\n"
             "3. 等差数列：[1,64,8]（起始,结束,步长，1-based）\n"
             "多个规范用空格分隔，整体用引号括起来"
    )
    parser.add_argument(
        '--splitSample', '-s', required=True,
        help='与splitCB数量对应的样本名称列表，顺序需严格匹配，用空格分隔，整体用引号括起来'
    )
    parser.add_argument(
        '--output', '-o', required=True,
        help='输出文件路径，带SampleID列的Cell_Summary_withSampleID.xls')
    
    args = parser.parse_args()
    
    # 手动分割参数为列表
    args.splitCB = args.splitCB.split()
    args.splitSample = args.splitSample.split()
    
    return args


def load_barcodes(path):
    """
    读取mapping.txt文件，提取barcode及其对应的0-based索引
    返回一个元组：(barcode_to_index, max_index)
        barcode_to_index: dict，key为barcode字符串，value为其对应的索引
        max_index: 最大索引值，用于验证范围有效性
    """
    barcode_to_index = {}
    # 正则表达式用于从第二列提取类似"bc2-xxx"中的数字部分
    pattern = re.compile(r'bc2-(\d+)')
    
    try:
        with open(path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue  # 跳过空行
                
                # 分割Tab分隔的列
                parts = line.split('\t')
                if len(parts) < 2:
                    print(f"WARNING: 第{line_num}行格式不正确，跳过该行", file=sys.stderr)
                    continue
                
                barcode = parts[0]
                info = parts[1]
                
                # 提取索引
                match = pattern.search(info)
                if not match:
                    print(f"WARNING: 第{line_num}行未找到有效的索引信息，跳过该行", file=sys.stderr)
                    continue
                
                try:
                    index = int(match.group(1))
                    barcode_to_index[barcode] = index
                except ValueError:
                    print(f"WARNING: 第{line_num}行索引不是有效整数，跳过该行", file=sys.stderr)
                    continue
    
    except IOError as e:
        print(f"ERROR: 无法读取文件{path}: {e}", file=sys.stderr)
        sys.exit(1)
    
    if not barcode_to_index:
        print(f"ERROR: mapping文件{path}中未找到有效的barcode和索引信息", file=sys.stderr)
        sys.exit(1)
    
    # 计算最大索引值
    max_index = max(barcode_to_index.values())
    return barcode_to_index, max_index



def parse_index_spec(spec, max_allowed_index):
    """
    解析单个索引规范，使用1-based索引（索引从1开始）
    支持三种格式：
    1. 范围格式：1-64 → 转换为[1,2,3,...,64]
    2. 单个索引：3 → 转换为[3]
    3. 等差数列：[1,64,8] → 1-based等差数列，包含起始和结束
    """
    indices = []

    # 处理等差数列格式 [start,end,step]
    if spec.startswith('[') and spec.endswith(']'):
        try:
            content = spec[1:-1].strip()
            # 分割并清理每个部分
            parts = [p.strip() for p in content.split(',')]
            if len(parts) != 3:
                raise ValueError("等差数列必须包含三个值: 起始,结束,步长")
                
            start, end, step = map(int, parts)
            # 校验步长
            if step <= 0:
                raise ValueError(f"步长必须为正整数，当前值：{step}")
            # 处理顺序问题，确保start <= end
            if start > end:
                start, end = end, start
            # 生成1-based索引
            current = start
            while current <= end:
                indices.append(current)
                current += step
        except Exception as e:
            print(f"ERROR: 无效的等差数列格式 '{spec}': {e}", file=sys.stderr)
            sys.exit(1)
    else:
        # 处理范围或单个索引格式
        parts = spec.split(',')
        for part in parts:
            part = part.strip()
            if '-' in part:
                # 范围格式
                try:
                    a, b = map(int, part.split('-'))
                    # 重要变更：使用1-based索引，范围检查从1开始
                    if a < 1 or b < a or b > max_allowed_index:
                        raise ValueError(f"范围超出有效索引范围(1-{max_allowed_index})")
                    # 生成1-based范围
                    indices.extend(range(a, b + 1))  # 包含结束值
                except ValueError as ve:
                    print(f"ERROR: 无效的范围格式 '{part}': {ve}", file=sys.stderr)
                    sys.exit(1)
                except:
                    print(f"ERROR: 无效的范围格式 '{part}'，应为1-based的start-end", file=sys.stderr)
                    sys.exit(1)
            else:
                # 单个索引
                try:
                    idx = int(part)
                    if idx < 1 or idx > max_allowed_index:
                        raise ValueError(f"索引超出有效范围(1-{max_allowed_index})")
                    indices.append(idx)
                except ValueError as ve:
                    print(f"ERROR: 无效的单个索引 '{part}': {ve}", file=sys.stderr)
                    sys.exit(1)
                except:
                    print(f"ERROR: 无效的单个索引 '{part}'，应为1-based整数", file=sys.stderr)
                    sys.exit(1)

    # 去重并保持顺序
    seen = set()
    unique_indices = [x for x in indices if x not in seen and not seen.add(x)]
    return unique_indices


def build_mapping(barcode_to_index, specs, samples, max_index):
    """
    根据specs和samples构建barcode->SampleID映射字典
    """
    # 首先创建索引到样本的映射
    index_to_sample = {}
    
    # 验证specs与samples长度一致
    if len(specs) != len(samples):
        print(f"ERROR: 提供的范围数量({len(specs)})与样本名称数量({len(samples)})不一致", file=sys.stderr)
        print(f"splitCB: {specs}")
        print(f"splitSample: {samples}")
        sys.exit(1)

    # 遍历每个索引规范，将对应索引映射到样本
    for spec, sample in zip(specs, samples):
        indices = parse_index_spec(spec, max_index)
        for idx in indices:
            if idx in index_to_sample:
                print(f"WARNING: 索引{idx}被分配到多个样本，使用最后一个分配的样本{sample}", file=sys.stderr)
            index_to_sample[idx] = sample
    
    # 然后创建barcode到样本的映射
    barcode_to_sample = {}
    for barcode, idx in barcode_to_index.items():
        barcode_to_sample[barcode] = index_to_sample.get(idx, 'Unknown')
    
    return barcode_to_sample


def assign_sample_id(df, mapping):
    """
    针对DataFrame的CB列，直接使用barcode匹配并映射SampleID
    """
    # 将CB列转换为字符串类型以确保正确匹配
    cb = df['CB'].astype(str)
    # 使用映射表进行转换，未知barcode标记为'Unknown'
    df['SampleID'] = cb.map(mapping).fillna('Unknown')
    return df


def main():
    # 解析命令行参数
    args = parse_args()

    # 读取barcode列表及其索引
    print(f"读取barcode信息: {args.barcodelist}")
    barcode_to_index, max_index = load_barcodes(args.barcodelist)
    print(f"成功加载{len(barcode_to_index)}个barcode，最大索引值为{max_index}")

    # 构建barcode到样本的映射表
    print("构建样本映射表...")
    print(f"splitCB参数: {args.splitCB}")
    print(f"splitSample参数: {args.splitSample}")
    mapping = build_mapping(barcode_to_index, args.splitCB, args.splitSample, max_index)

    # 读取输入表格
    try:
        print(f"读取输入文件: {args.inputfile}")
        df = pd.read_csv(args.inputfile, sep='\t', dtype=str)
    except Exception as e:
        print(f"ERROR: 读取文件失败 {args.inputfile}: {e}", file=sys.stderr)
        sys.exit(1)

    # 检查CB列是否存在
    if 'CB' not in df.columns:
        print("ERROR: 输入文件中未检测到CB列", file=sys.stderr)
        sys.exit(1)

    # 添加SampleID列
    print("分配样本ID...")
    df = assign_sample_id(df, mapping)

    # 统计各样本数量
    sample_counts = df['SampleID'].value_counts()
    print("样本分配统计:")
    for sample, count in sample_counts.items():
        print(f"  {sample}: {count}个细胞")

    # 写出结果
    try:
        df.to_csv(args.output, sep='\t', index=False)
        print(f"成功: 输出文件已保存至 {args.output}")
    except Exception as e:
        print(f"ERROR: 写入文件失败 {args.output}: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
