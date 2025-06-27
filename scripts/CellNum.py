# -*- coding: utf-8 -*-
"""
CellNum.py (Fixed)

主要：
1. 修正等差数列语法解析逻辑，现在可以正确识别 [start,end,step] 格式
2. 优化参数处理顺序，避免逗号分割干扰等差数列解析
"""

import argparse
import gzip
import os
import sys
from pathlib import Path
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser(description="按 CB2 序列拆分混样（修复版）")
    p.add_argument(
        "--inputfileCellMatrix", 
        nargs="+", required=True,
        help="所有 filtered_feature_bc_matrix 文件夹路径列表"
    )
    p.add_argument(
        "--inputfileBarcodelist", 
        required=True,
        help="CB2 序列列表文件"
    )
    p.add_argument(
        "--splitCB",
        nargs="+", required=True,
        help="索引规范，支持：1-5, 1,3,5, [1,89,8]"
    )
    p.add_argument(
        "--splitSample",
        nargs="+", required=True,
        help="样本名称列表"
    )
    p.add_argument(
        "--outCellNum",
        required=True,
        help="输出统计表格路径"
    )
    p.add_argument(
        "--outCellList",
        nargs="+", required=True,
        help="输出barcode列表文件"
    )
    return p.parse_args()

def load_barcode_list(barcode_list_file):
    """读取并返回barcode列表"""
    with open(barcode_list_file) as f:
        return [line.strip() for line in f if line.strip()]

def parse_index_spec(spec):
    """
    完全重构的解析逻辑：
    1. 优先处理等差数列语法
    2. 再处理常规语法
    """
    indices = []
    
    # ===== 新增：处理完整等差数列语法 =====
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
    
    # ===== 原始解析逻辑 =====
    parts = [p.strip() for p in spec.split(',') if p.strip()]
    
    for part in parts:
        if '-' in part:
            try:
                a, b = map(int, part.split('-'))
                a, b = sorted([a, b])
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

def main():
    args = parse_args()
    
    # 参数校验
    if len(args.splitCB) != len(args.splitSample):
        sys.exit("错误: --splitCB 和 --splitSample 数量不匹配")
    if len(args.outCellList) != len(args.splitSample):
        sys.exit("错误: --outCellList 数量不匹配")

    # 创建输出目录（使用Pathlib替代os.makedirs）
    output_dir = Path("06-CB2_Sample")
    output_dir.mkdir(exist_ok=True)
    
    # 加载barcode列表
    all_barcodes = load_barcode_list(args.inputfileBarcodelist)
    total = len(all_barcodes)
    if total == 0:
        sys.exit("错误: Barcode列表为空")

    # 解析索引规范
    sample_indices = {}
    for spec, sample in zip(args.splitCB, args.splitSample):
        try:
            indices = parse_index_spec(spec)
            # 校验索引范围
            for idx in indices:
                if idx < 0 or idx >= total:
                    sys.exit(f"错误: 索引 {idx+1} 超出范围 (1-{total})")
            sample_indices[sample] = indices
        except ValueError as e:
            sys.exit(f"参数解析错误: {e}")

    # 准备输出文件
    sample_files = {
        sample: open(output_dir / path, 'w')  # 使用Pathlib拼接路径
        for sample, path in zip(args.splitSample, args.outCellList)
    }
    counts = {sample:0 for sample in args.splitSample}

    # 处理矩阵目录
    for matrix_dir in args.inputfileCellMatrix:
        gz_path = os.path.join(matrix_dir, "barcodes.tsv.gz")
        if not os.path.exists(gz_path):
            sys.exit(f"错误: 文件不存在 {gz_path}")
            
        with gzip.open(gz_path, 'rt') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                    
                # 提取第二段barcode
                parts = line.split('_')
                if len(parts) < 2:
                    continue
                cb2 = parts[1]
                
                # 匹配样本
                for sample, indices in sample_indices.items():
                    if cb2 in {all_barcodes[i] for i in indices}:
                        counts[sample] += 1
                        sample_files[sample].write(line + '\n')
                        break

    # 关闭所有文件
    for f in sample_files.values():
        f.close()

    # 输出结果
    df = pd.DataFrame.from_dict(counts, orient='index', columns=['Estimated Cells'])
    dir_path = '06-CB2_Sample'
    os.makedirs(dir_path, exist_ok=True)
    full_path = os.path.join(dir_path, args.outCellNum)
    df.T.to_csv(full_path, sep='\t', index=False)
    
    print(">> 处理完成：")
    print(df.to_string())

if __name__ == "__main__":
    main()