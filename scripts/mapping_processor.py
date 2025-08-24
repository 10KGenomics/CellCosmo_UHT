#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mapping_processor.py (Optimized)

功能：
1. 从mapping.txt文件中读取barcode和bc2-*信息
2. 支持多种分类规则：区间、等差数列、逗号分隔列表
3. 生成样本对应的barcode列表文件
4. 统计EstimatedCell与AllCell的匹配百分比

分类规则支持：
- 区间: 1-64
- 等差数列: [1,89,8]
- 逗号分隔列表: 1,2,3,4
- 混合模式: 1-5,7,9,11-13
"""

import re
import os
import sys
import argparse
from collections import defaultdict

def parse_index_spec(spec):
    """
    解析索引规范：
    支持区间、等差数列和逗号分隔列表
    返回0-based索引列表
    """
    indices = []
    
    # 处理等差数列语法 [start,end,step]
    if spec.startswith('[') and spec.endswith(']'):
        try:
            content = spec[1:-1].strip()
            params = [p.strip() for p in content.split(',')]
            if len(params) != 3:
                raise ValueError("需要 exactly 3 个参数")
                
            start = int(params[0])
            end = int(params[1])
            step = int(params[2])
            
            if step <= 0:
                raise ValueError("步长必须 >0")
                
            # 自动调整方向
            if start > end:
                start, end = end, start
                step = abs(step)
            elif start < end and step < 0:
                step = abs(step)
                
            current = start
            while current <= end if start <= end else current >= end:
                indices.append(current - 1)  # 转换为 0-based
                current += step
                
            return list(dict.fromkeys(indices))  # 去重保序
            
        except Exception as e:
            raise ValueError(f"无效的等差数列语法 '{spec}': {e}")
    
    # 处理常规语法 (x-y 或 x,y,z)
    parts = [p.strip() for p in spec.split(',') if p.strip()]
    
    for part in parts:
        if '-' in part:
            try:
                a, b = map(int, part.split('-'))
                a, b = sorted([a, b])
                # 扩展范围：a-1 到 b（包含b）
                indices.extend(range(a-1, b))
            except:
                raise ValueError(f"无效区间格式 '{part}'")
        else:
            try:
                indices.append(int(part)-1)
            except:
                raise ValueError(f"无效索引值 '{part}'")
    
    # 去重并保持顺序
    seen = set()
    return [x for x in indices if not (x in seen or seen.add(x))]

def parse_mapping_file(mapping_path):
    """
    解析mapping.txt文件
    返回: 
      barcode_data: 列表，包含所有有效的(barcode, bc2_value)元组
      error_count: 解析错误的数量
    """
    barcode_data = []
    error_count = 0
    bc2_pattern = re.compile(r'bc2-(\d+)')
    
    try:
        with open(mapping_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split('\t')
                if len(parts) < 2:
                    error_count += 1
                    continue
                
                barcode = parts[0].strip()
                info_str = parts[1].strip()
                
                # 从信息字符串中提取bc2值
                match = bc2_pattern.search(info_str)
                if not match:
                    error_count += 1
                    continue
                
                try:
                    bc2_value = int(match.group(1))
                    barcode_data.append((barcode, bc2_value))
                except ValueError:
                    error_count += 1
                    
    except FileNotFoundError:
        print(f"错误: 文件 '{mapping_path}' 不存在")
        sys.exit(1)
    
    return barcode_data, error_count

def compare_cell_lists(all_cell_files, estimated_cell_dir, prefix="3-"):
    """
    比较AllCell和EstimatedCell列表
    返回每个样本的匹配百分比
    """
    comparison_results = {}
    
    for sample, all_cell_file in all_cell_files.items():
        # 读取AllCell列表
        try:
            with open(all_cell_file, 'r') as f:
                all_cell_barcodes = set(line.strip() for line in f)
        except FileNotFoundError:
            print(f"警告: AllCell文件 {all_cell_file} 不存在")
            continue
            
        # 读取EstimatedCell列表
        estimated_file = os.path.join(estimated_cell_dir, f"{prefix}EstimatedCell_{sample}.list")
        try:
            with open(estimated_file, 'r') as f:
                estimated_barcodes = set(line.strip() for line in f)
        except FileNotFoundError:
            print(f"警告: EstimatedCell文件 {estimated_file} 不存在")
            continue
            
        # 计算匹配的barcode
        matched_barcodes = all_cell_barcodes & estimated_barcodes
        total_estimated = len(estimated_barcodes)
        match_count = len(matched_barcodes)
        
        # 计算匹配百分比
        match_percent = (match_count / total_estimated * 100) if total_estimated > 0 else 0
        
        comparison_results[sample] = {
            'total_estimated': total_estimated,
            'match_count': match_count,
            'match_percent': match_percent
        }
    
    return comparison_results

def main():
    # 设置命令行参数
    parser = argparse.ArgumentParser(description='Barcode分类与统计工具')
    parser.add_argument('--input', required=True, help='mapping.txt文件路径')
    parser.add_argument('--splitCB', nargs='+', required=True, help='分类规则列表')
    parser.add_argument('--samples', nargs='+', required=True, help='样本名称列表')
    parser.add_argument('--output', nargs='+', required=True, help='输出文件列表')
    parser.add_argument('--estimated-dir', default='06-CB2_Sample', 
                        help='EstimatedCell文件目录 (default: 06-CB2_Sample)')
    parser.add_argument('--prefix', default='3-', 
                        help='文件前缀 (default: 3-)')
    args = parser.parse_args()
    
    # 校验参数
    if len(args.splitCB) != len(args.samples):
        print("错误: --splitCB 和 --samples 数量不匹配")
        sys.exit(1)
    if len(args.output) != len(args.samples):
        print("错误: --output 数量不匹配")
        sys.exit(1)
    
    # 解析分类规则
    sample_rules = {}
    sample_indices = {}
    print("\n>> 解析分类规则:")
    for i, (spec, sample) in enumerate(zip(args.splitCB, args.samples)):
        try:
            indices = parse_index_spec(spec)
            sample_rules[sample] = spec
            sample_indices[sample] = set(indices)  # 转换为集合便于快速查找
            print(f"  样本 {sample}: 规则 '{spec}' -> {len(indices)} 个索引")
        except ValueError as e:
            print(f"错误: 解析规则 '{spec}' 失败 - {e}")
            sys.exit(1)
    
    # 加载mapping数据
    print("\n>> 加载mapping文件...")
    barcode_data, error_count = parse_mapping_file(args.input)
    total_entries = len(barcode_data) + error_count
    
    if not barcode_data:
        print("错误: 没有找到有效的barcode数据")
        sys.exit(1)
    
    print(f"  处理完成: 共解析 {len(barcode_data)} 个有效条目")
    if error_count > 0:
        print(f"  警告: 跳过 {error_count} 个无效条目 ({error_count/total_entries*100:.2f}%)")
    
    # 分类barcode到样本
    print("\n>> 按样本分类barcode...")
    sample_barcodes = defaultdict(list)
    out_of_range_count = 0
    
    for barcode, bc2_value in barcode_data:
        # 转换为0-based索引
        bc2_index0 = bc2_value - 1
        
        # 分配到样本
        assigned = False
        for sample, indices_set in sample_indices.items():
            if bc2_index0 in indices_set:
                sample_barcodes[sample].append(barcode)
                assigned = True
                break
                
        if not assigned:
            out_of_range_count += 1
    
    # 输出统计信息
    print("\n>> 分类统计:")
    total_classified = 0
    for sample in args.samples:
        count = len(sample_barcodes.get(sample, []))
        total_classified += count
        print(f"  样本 {sample}: {count} 个barcode")
    
    print(f"\n  总计分类: {total_classified} 个barcode")
    if out_of_range_count > 0:
        print(f"  警告: {out_of_range_count} 个barcode的bc2值不在任何样本范围内")
    
    # 写入输出文件
    print("\n>> 写入输出文件...")
    all_cell_files = {}
    for sample, output_file in zip(args.samples, args.output):
        barcodes = sample_barcodes.get(sample, [])
        try:
            with open(output_file, 'w') as f:
                for barcode in barcodes:
                    f.write(barcode + '\n')
            print(f"  已写入 {len(barcodes)} 个barcode到 {output_file}")
            all_cell_files[sample] = output_file
        except IOError as e:
            print(f"错误: 写入文件 {output_file} 失败 - {e}")
    
    # 比较AllCell和EstimatedCell列表
    print("\n>> 比较AllCell和EstimatedCell列表...")
    comparison = compare_cell_lists(
        all_cell_files, 
        args.estimated_dir,
        prefix=args.prefix
    )
    
    # 输出比较结果
    print("\n>> 匹配统计:")
    print(f"{'样本':<8}{'Estimated Cells':>15}{'Matched Cells':>15}{'Match %':>10}")
    total_estimated = 0
    total_matched = 0
    
    for sample in args.samples:
        if sample in comparison:
            stats = comparison[sample]
            total_estimated += stats['total_estimated']
            total_matched += stats['match_count']
            print(f"{sample:<8}{stats['total_estimated']:>15}{stats['match_count']:>15}{stats['match_percent']:>10.2f}%")
        else:
            print(f"{sample:<8}{'N/A':>15}{'N/A':>15}{'N/A':>10}")
    
    # 计算总体匹配率
    if total_estimated > 0:
        overall_percent = total_matched / total_estimated * 100
        print(f"\n  总体匹配率: {overall_percent:.2f}% ({total_matched}/{total_estimated})")
    else:
        print("\n  警告: 没有可用的EstimatedCell数据")
    
    print("\n>> 处理完成!")

if __name__ == "__main__":
    main()