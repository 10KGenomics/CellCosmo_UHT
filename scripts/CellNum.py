# -*- coding: utf-8 -*-
"""
CellNum.py (Enhanced)

主要优化：
1. 支持新的 mapping.txt 格式（包含16bp barcode和文库信息）
2. 解析 bc2-xxx 格式提取CB2索引号
3. 优化内存使用，避免加载整个mapping.txt到内存
4. 增加详细的错误处理和日志
5. 改进索引分配逻辑，支持任意索引范围
6. 添加类型注解和详细注释
"""

import argparse
import gzip
import os
import re
import sys
from collections import defaultdict
import pandas as pd
import logging
from typing import Dict, List, Set, Tuple

# 配置日志系统
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger(__name__)

def parse_args() -> argparse.Namespace:
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description="按 CB2 索引拆分混样（增强版）",
        epilog="示例: python CellNum.py --inputfileCellMatrix matrix_dir1 matrix_dir2 "
               "--inputfileBarcodelist mapping.txt --splitCB '1-64 65-128 129-192' "
               "--splitSample A B C --outCellNum cell_counts.tsv "
               "--outCellList cellA.list cellB.list cellC.list"
    )
    parser.add_argument(
        "--inputfileCellMatrix", 
        nargs="+", required=True,
        help="所有 filtered_feature_bc_matrix 文件夹路径列表"
    )
    parser.add_argument(
        "--inputfileBarcodelist", 
        required=True,
        help="mapping.txt 文件路径（包含barcode和文库信息）"
    )
    parser.add_argument(
        "--splitCB",
        nargs="+", required=True,
        help="CB2索引范围规范，如: '1-64 65-128 129-192'"
    )
    parser.add_argument(
        "--splitSample",
        nargs="+", required=True,
        help="样本名称列表，如: 'A B C'"
    )
    parser.add_argument(
        "--outCellNum",
        required=True,
        help="输出细胞计数表格路径"
    )
    parser.add_argument(
        "--outCellList",
        nargs="+", required=True,
        help="输出barcode列表文件路径列表"
    )
    return parser.parse_args()

def parse_index_spec(spec: str) -> Set[int]:
    """
    解析索引范围规范，返回索引集合
    
    支持格式:
    - 单个数字: "1"
    - 范围: "1-64"
    - 混合: "1,3,5-10"
    - 等差数列: "[1,100,5]"
    """
    indices = set()
    
    # 处理等差数列语法 [start,end,step]
    if spec.startswith('[') and spec.endswith(']'):
        try:
            parts = spec[1:-1].split(',')
            if len(parts) != 3:
                raise ValueError("需要3个参数: start,end,step")
                
            start = int(parts[0].strip())
            end = int(parts[1].strip())
            step = int(parts[2].strip())
            
            if step == 0:
                raise ValueError("步长不能为0")
                
            # 生成数列
            current = start
            if start <= end:
                while current <= end:
                    indices.add(current)
                    current += abs(step)
            else:
                while current >= end:
                    indices.add(current)
                    current -= abs(step)
                    
            return indices
            
        except Exception as e:
            raise ValueError(f"无效的等差数列语法 '{spec}': {str(e)}")
    
    # 处理常规语法 (逗号分隔的单个数字或范围)
    parts = [p.strip() for p in spec.split(',') if p.strip()]
    
    for part in parts:
        if '-' in part:
            try:
                a, b = map(int, part.split('-'))
                low, high = sorted([a, b])
                indices.update(range(low, high + 1))
            except:
                raise ValueError(f"无效范围格式 '{part}'")
        else:
            try:
                indices.add(int(part))
            except:
                raise ValueError(f"无效索引值 '{part}'")
    
    return indices

def parse_bc2_info(info_str: str) -> int:
    """
    从文库信息字符串中提取CB2索引号
    
    示例输入: "bc20,a5,l1,bc2-135"
    返回: 135
    """
    # 查找 bc2-xxx 模式
    match = re.search(r'bc2-(\d+)', info_str)
    if match:
        return int(match.group(1))
    
    # 查找其他可能格式
    match = re.search(r'bc2_(\d+)', info_str)
    if match:
        return int(match.group(1))
    
    raise ValueError(f"未找到有效的CB2索引号: {info_str}")

def create_sample_index_map(
    mapping_file: str, 
    split_specs: List[str], 
    samples: List[str]
) -> Tuple[Dict[str, Set[int]], Dict[str, str]]:
    """
    创建样本索引映射
    
    返回:
    - sample_indices: 每个样本的CB2索引集合 {样本名: {索引1, 索引2, ...}}
    - barcode_to_sample: barcode到样本的映射 {16bp_barcode: 样本名}
    """
    # 解析每个样本的索引范围
    if len(split_specs) != len(samples):
        raise ValueError("splitCB和splitSample数量不匹配")
    
    sample_indices = {}
    for spec, sample in zip(split_specs, samples):
        try:
            sample_indices[sample] = parse_index_spec(spec)
        except ValueError as e:
            raise ValueError(f"样本 '{sample}' 的索引规范错误: {str(e)}")
    
    # 创建barcode到样本的映射
    barcode_to_sample = {}
    bc2_index_count = defaultdict(int)
    
    try:
        with open(mapping_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    logger.warning(f"跳过无效行 {line_num}: {line.strip()}")
                    continue
                
                barcode = parts[0].strip()
                info_str = parts[1].strip()
                
                try:
                    # 解析CB2索引号
                    bc2_index = parse_bc2_info(info_str)
                    bc2_index_count[bc2_index] += 1
                    
                    # 查找对应的样本
                    for sample, indices in sample_indices.items():
                        if bc2_index in indices:
                            barcode_to_sample[barcode] = sample
                            break
                    else:
                        # 如果没找到匹配的样本
                        logger.debug(f"CB2索引 {bc2_index} 不属于任何样本范围")
                
                except ValueError as e:
                    logger.warning(f"行 {line_num} 解析错误: {str(e)}")
    
    except IOError as e:
        raise IOError(f"无法读取mapping文件: {str(e)}")
    
    # 检查是否有barcode被分配
    if not barcode_to_sample:
        raise RuntimeError("没有barcode被分配到样本，请检查mapping.txt和splitCB参数")
    
    # 记录统计信息
    logger.info(f"成功映射 {len(barcode_to_sample)} 个barcode到样本")
    logger.debug(f"CB2索引分布: {dict(bc2_index_count)}")
    
    return sample_indices, barcode_to_sample

def main():
    args = parse_args()
    
    # 参数验证
    if len(args.splitSample) != len(args.outCellList):
        logger.error("--splitSample 和 --outCellList 数量不匹配")
        sys.exit(1)
    
    # 创建输出目录
    output_dir = "06-CB2_Sample"
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        # 步骤1: 创建样本索引映射
        sample_indices, barcode_to_sample = create_sample_index_map(
            args.inputfileBarcodelist,
            args.splitCB,
            args.splitSample
        )
    except Exception as e:
        logger.error(f"映射创建失败: {str(e)}")
        sys.exit(1)
    
    # 步骤2: 准备输出文件和计数器
    sample_files = {}
    sample_counts = {sample: 0 for sample in args.splitSample}
    
    for sample, out_path in zip(args.splitSample, args.outCellList):
        full_path = os.path.join(output_dir, out_path)
        sample_files[sample] = open(full_path, 'w')
    
    # 步骤3: 处理所有矩阵目录
    total_barcodes = 0
    matched_barcodes = 0
    
    for matrix_dir in args.inputfileCellMatrix:
        # 查找barcodes文件 (支持 .gz 和未压缩)
        barcodes_path = None
        for fname in ["barcodes.tsv", "barcodes.tsv.gz"]:
            path = os.path.join(matrix_dir, fname)
            if os.path.exists(path):
                barcodes_path = path
                break
        
        if not barcodes_path:
            logger.error(f"在 {matrix_dir} 中找不到barcodes文件")
            continue
        
        logger.info(f"处理矩阵目录: {matrix_dir}")
        
        # 打开barcodes文件 (支持gzip)
        open_func = gzip.open if barcodes_path.endswith('.gz') else open
        mode = 'rt' if barcodes_path.endswith('.gz') else 'r'
        
        try:
            with open_func(barcodes_path, mode) as f:
                for line in f:
                    total_barcodes += 1
                    barcode = line.strip()
                    
                    # 在mapping中查找样本
                    sample = barcode_to_sample.get(barcode)
                    if sample:
                        matched_barcodes += 1
                        sample_counts[sample] += 1
                        sample_files[sample].write(line)
        except Exception as e:
            logger.error(f"处理 {barcodes_path} 失败: {str(e)}")
    
    # 步骤4: 关闭所有输出文件
    for f in sample_files.values():
        f.close()
    
    # 步骤5: 输出统计结果
    logger.info(f"处理完成: 共扫描 {total_barcodes} 个barcode, 匹配 {matched_barcodes} 个")
    
    # 创建计数DataFrame
    df = pd.DataFrame.from_dict(
        sample_counts, 
        orient='index', 
        columns=['Cells Num']
    )
    
    # 保存到文件
    output_path = os.path.join(output_dir, args.outCellNum)
    df.to_csv(output_path, sep='\t', index_label="Sample")
    
    # 打印结果
    logger.info(f"细胞计数已保存至: {output_path}")
    logger.info("\n" + df.to_string())
    
    # 记录输出文件路径
    for sample, out_path in zip(args.splitSample, args.outCellList):
        logger.info(f"样本 {sample} barcode列表: {os.path.join(output_dir, out_path)}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logger.exception("处理过程中发生致命错误")
        sys.exit(1)