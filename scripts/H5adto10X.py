# -*- coding: utf-8 -*-
"""
h5adto10X.py
-------------
该脚本用于将 h5ad 格式的单细胞数据转换为 10X Genomics 标准格式，即生成包含
  - barcodes.tsv.gz
  - features.tsv.gz
  - matrix.mtx.gz
三个文件的目录（例如：SampleA_filtered_feature_bc_matrix）
用法:
    python h5adto10X.py SampleA SampleA_filtered_feature_bc_matrix
其中 SampleA.h5ad 为输入文件，转换结果保存在指定的输出目录内。
"""

import os
import sys
import gzip
import shutil
import scanpy as sc
import numpy as np
from scipy.io import mmwrite
from scipy import sparse

def usage():
    print("Usage: python h5adto10X.py <sample_basename> <output_directory>")
    print("Example: python h5adto10X.py SampleA SampleA_filtered_feature_bc_matrix")
    sys.exit(1)

def export_barcodes(adata, out_path):
    """
    将细胞条形码（barcode）写入压缩后的 barcodes.tsv.gz 文件。
    每行包含一个细胞的 barcode，来源于 AnnData 对象的 obs_names。
    """
    # 使用 gzip.open 写入文本信息
    with gzip.open(out_path, 'wt') as f:
        for barcode in adata.obs_names:
            f.write(f"{barcode}\n")
    print(f"Barcodes written to {out_path}")

def export_features(adata, out_path):
    """
    将基因（features）信息写入压缩后的 features.tsv.gz 文件。
    输出文件包含三列：
      1. gene_id（采用 adata.var_names）
      2. gene_name（如果 adata.var 中有 'gene_symbols' 列则使用之，否则使用 var_names）
      3. feature_type，统一设定为 "Gene Expression"
    """
    # 检查 adata.var 中是否包含 gene_symbols 字段
    if 'gene_symbols' in adata.var.columns:
        gene_names = adata.var['gene_symbols'].tolist()
    else:
        gene_names = adata.var_names.tolist()
    
    with gzip.open(out_path, 'wt') as f:
        # 对于每个基因写出三列内容，tab 分隔
        for gene_id, gene_name in zip(adata.var_names, gene_names):
            f.write(f"{gene_id}\t{gene_name}\tGene Expression\n")
    print(f"Features written to {out_path}")

def export_matrix(adata, out_path):
    """
    将 AnnData 对象中的表达矩阵转换为 Matrix Market 格式，写入 matrix.mtx.gz 文件。
    注意：10X Genomics 格式的 matrix.mtx 中行代表基因、列代表细胞，
    因此需要对 AnnData 中 (cells, genes) 的矩阵进行转置。
    """
    # 由于 mmwrite 无法直接写入 gzip 文件，先写入临时文件，再压缩后删除临时文件
    temp_matrix = out_path.replace(".gz", "")
    
    # 如果表达矩阵为稠密矩阵，也转换为稀疏矩阵
    if not sparse.isspmatrix(adata.X):
        X = sparse.coo_matrix(adata.X)
    else:
        # 转置矩阵，使其形状为 (n_genes, n_cells)
        X = adata.X.transpose().tocoo()
    
    # 使用 mmwrite 写入临时文件
    mmwrite(temp_matrix, X)
    print(f"Matrix Market file written to temporary file {temp_matrix}")
    
    # 将临时文件压缩
    with open(temp_matrix, 'rb') as f_in, gzip.open(out_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    print(f"Compressed matrix written to {out_path}")
    
    # 删除临时文件
    os.remove(temp_matrix)

def main():
    # 判断命令行参数数量
    if len(sys.argv) != 3:
        usage()
    
    sample_name = sys.argv[1]         # 例如 "SampleA"
    output_dir = sys.argv[2]          # 例如 "SampleA_filtered_feature_bc_matrix"
    
    # 构建输入 h5ad 文件路径，默认文件名为 SampleA.h5ad
    input_h5ad = f"06-CB2_Sample/{sample_name}.h5ad"
    if not os.path.exists(input_h5ad):
        print(f"Error: Cannot find {input_h5ad}")
        sys.exit(1)
    
    print(f"Loading {input_h5ad} ...")
    adata = sc.read(input_h5ad)
    print("h5ad file loaded successfully.")
    
    # 确保输出目录存在，若不存在则创建之
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created directory {output_dir}")
    
    # 定义输出文件路径
    barcodes_out = os.path.join(output_dir, "barcodes.tsv.gz")
    features_out = os.path.join(output_dir, "features.tsv.gz")
    matrix_out   = os.path.join(output_dir, "matrix.mtx.gz")
    
    # 输出 10X Genomics 的三个文件
    export_barcodes(adata, barcodes_out)
    export_features(adata, features_out)
    export_matrix(adata, matrix_out)
    
    print("Conversion completed successfully.")

if __name__ == "__main__":
    main()
