# -*- coding: utf-8 -*-
"""
h5adto10X.py
-------------
优化版脚本：将 h5ad 格式的单细胞数据转换为 10X Genomics 标准格式

主要改进：
1. 输入文件路径作为命令行参数（不再硬编码）
2. 输出目录作为命令行参数
3. 增加参数校验和错误处理
4. 添加详细的进度日志
5. 支持多种稀疏矩阵格式

用法:
    python h5adto10X.py <input_h5ad> <output_directory>
    
示例:
    # 处理过滤矩阵
    python h5adto10X.py \\
        "06-CB2_Sample/filtered/A_filtered.h5ad" \\
        "07-CB2_Matrix/A_filtered_feature_bc_matrix"
    
    # 处理原始矩阵
    python h5adto10X.py \\
        "06-CB2_Sample/raw/A_raw.h5ad" \\
        "07-CB2_Matrix/A_raw_feature_bc_matrix"
"""

import os
import sys
import gzip
import shutil
import scanpy as sc
import numpy as np
from scipy.io import mmwrite
from scipy import sparse
import argparse

def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='将 h5ad 格式转换为 10X Genomics 标准格式',
        epilog='示例: python h5adto10X.py input.h5ad output_directory'
    )
    parser.add_argument(
        'input_h5ad', 
        type=str,
        help='输入h5ad文件路径 (e.g., "06-CB2_Sample/filtered/A_filtered.h5ad")'
    )
    parser.add_argument(
        'output_dir', 
        type=str,
        help='输出目录路径 (e.g., "07-CB2_Matrix/A_filtered_feature_bc_matrix")'
    )
    return parser.parse_args()

def validate_paths(input_h5ad, output_dir):
    """验证输入文件存在且输出目录可创建"""
    # 检查输入文件
    if not os.path.isfile(input_h5ad):
        raise FileNotFoundError(f"输入文件不存在: {input_h5ad}")
    
    # 检查文件扩展名
    if not input_h5ad.lower().endswith('.h5ad'):
        print(f"警告: 输入文件 '{os.path.basename(input_h5ad)}' 可能不是h5ad格式")
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    print(f"输出目录已创建: {output_dir}")

def export_barcodes(adata, out_path):
    """
    将细胞条形码（barcode）写入压缩后的 barcodes.tsv.gz 文件
    - 每行包含一个细胞的 barcode
    - 来源: AnnData 对象的 obs_names
    """
    # 获取barcode列表
    barcodes = adata.obs_names.tolist()
    
    # 写入gzip压缩文件
    with gzip.open(out_path, 'wt') as f:
        for barcode in barcodes:
            f.write(f"{barcode}\n")
    print(f"✓ 成功写入 {len(barcodes)} 个barcodes到 {out_path}")

def export_features(adata, out_path):
    """
    将基因信息写入压缩后的 features.tsv.gz 文件
    格式:
      gene_id <tab> gene_symbol <tab> feature_type
    
    字段优先级:
      1. gene_id: 使用adata.var_names
      2. gene_symbol: 
          - 优先使用adata.var中的'gene_symbols'列
          - 其次使用'gene_names'列
          - 最后使用gene_id
      3. feature_type: 统一为"Gene Expression"
    """
    # 确定基因符号列
    symbol_col = None
    for col in ['gene_symbols', 'gene_names', 'gene_symbol', 'symbol']:
        if col in adata.var.columns:
            symbol_col = col
            break
    
    # 使用可用的列或回退到gene_id
    if symbol_col:
        gene_symbols = adata.var[symbol_col].tolist()
        print(f"  使用 '{symbol_col}' 列作为基因符号")
    else:
        gene_symbols = adata.var_names.tolist()
        print("  警告: 未找到基因符号列，使用gene_id作为符号")
    
    # 写入文件
    with gzip.open(out_path, 'wt') as f:
        for gene_id, gene_name in zip(adata.var_names, gene_symbols):
            f.write(f"{gene_id}\t{gene_name}\tGene Expression\n")
    
    print(f"✓ 成功写入 {len(adata.var_names)} 个features到 {out_path}")

def export_matrix(adata, out_path):
    """
    将表达矩阵转换为 Matrix Market 格式并压缩
    
    处理流程:
      1. 转置矩阵: (cells × genes) → (genes × cells)
      2. 转换为COO稀疏格式
      3. 写入临时.mtx文件
      4. 压缩为.mtx.gz
      5. 删除临时文件
    
    支持矩阵类型:
      - CSR, CSC, COO, BSR (所有scipy稀疏格式)
      - 稠密矩阵(numpy.ndarray)
    """
    # 获取矩阵并转置
    X = adata.X
    
    # 处理稠密矩阵
    if isinstance(X, np.ndarray):
        print("  检测到稠密矩阵，转换为稀疏格式...")
        X = sparse.coo_matrix(X)
    
    # 确保是稀疏矩阵并转置
    if sparse.issparse(X):
        X = X.transpose().tocoo()
    else:
        raise TypeError(f"不支持的数据类型: {type(X)}")
    
    print(f"  矩阵形状: 基因={X.shape[0]}, 细胞={X.shape[1]}, 非零元素={X.nnz}")
    
    # 创建临时文件路径
    temp_matrix = os.path.join(os.path.dirname(out_path), "temp_matrix.mtx")
    
    # 写入Matrix Market格式
    mmwrite(temp_matrix, X)
    print(f"  临时矩阵文件已写入: {temp_matrix}")
    
    # 压缩文件
    with open(temp_matrix, 'rb') as f_in, gzip.open(out_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    
    # 删除临时文件
    os.remove(temp_matrix)
    print(f"✓ 成功写入并压缩矩阵到 {out_path}")

def main():
    """主函数：执行转换流程"""
    try:
        # 解析参数
        args = parse_arguments()
        input_h5ad = os.path.abspath(args.input_h5ad)
        output_dir = os.path.abspath(args.output_dir)
        
        print("=" * 60)
        print(f"开始转换: {os.path.basename(input_h5ad)} → 10X格式")
        print(f"输入文件: {input_h5ad}")
        print(f"输出目录: {output_dir}")
        print("=" * 60)
        
        # 验证路径
        validate_paths(input_h5ad, output_dir)
        
        # 加载h5ad文件
        print(f"加载h5ad文件: {input_h5ad}")
        adata = sc.read(input_h5ad)
        print(f"✓ 加载成功! 细胞数={adata.n_obs}, 基因数={adata.n_vars}")
        
        # 定义输出文件路径
        barcodes_out = os.path.join(output_dir, "barcodes.tsv.gz")
        features_out = os.path.join(output_dir, "features.tsv.gz")
        matrix_out   = os.path.join(output_dir, "matrix.mtx.gz")
        
        # 导出10X格式文件
        export_barcodes(adata, barcodes_out)
        export_features(adata, features_out)
        export_matrix(adata, matrix_out)
        
        # 生成完成报告
        print("\n" + "=" * 60)
        print("转换完成! 生成以下文件:")
        print(f"  - barcodes.tsv.gz: {os.path.getsize(barcodes_out)/1024/1024:.2f} MB")
        print(f"  - features.tsv.gz: {os.path.getsize(features_out)/1024/1024:.2f} MB")
        print(f"  - matrix.mtx.gz:   {os.path.getsize(matrix_out)/1024/1024:.2f} MB")
        print("=" * 60)
        
    except Exception as e:
        print(f"\n❌ 错误: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()