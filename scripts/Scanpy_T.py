import os
from glob import glob
import argparse
import scanpy as sc
import pandas as pd
import numpy as np

##########################################
# 定义一个函数用于计算 AnnData 对象的基本统计指标
##########################################
def compute_library_stats(adata):
    """
    计算单个 AnnData 对象的统计指标：
      - Num_Cells: 细胞数
      - Total_UMI: 全部转录本（UMI）数
      - Median_UMI_per_Cell: 每个细胞的 UMI 中位数
      - Mean_UMI_per_Cell: 每个细胞的 UMI 平均数
      - Median_Genes_per_Cell: 每个细胞检测到的基因数中位值
      - Mean_Genes_per_Cell: 每个细胞检测到的基因数均值
      - Total_Genes_Covered: 文库中总共检测到的基因种类数（至少在1个细胞中有表达）
    """
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    num_cells = adata.n_obs
    total_umi = adata.obs['total_counts'].sum()
    median_umi = np.median(adata.obs['total_counts'])
    mean_umi = np.mean(adata.obs['total_counts'])
    median_genes = np.median(adata.obs['n_genes_by_counts'])
    mean_genes = np.mean(adata.obs['n_genes_by_counts'])

    if hasattr(adata.X, 'toarray'):
        gene_totals = np.array(adata.X.toarray().sum(axis=0)).flatten()
    else:
        gene_totals = np.array(adata.X.sum(axis=0)).flatten()
    gene_covered = np.sum(gene_totals > 0)
    return {
        'Num_Cells': num_cells,
        'Total_UMI': total_umi,
        'Median_UMI_per_Cell': median_umi,
        'Mean_UMI_per_Cell': mean_umi,
        'Median_Genes_per_Cell': median_genes,
        'Mean_Genes_per_Cell': mean_genes,
        'Total_Genes_Covered': gene_covered
    }

def parse_args():
    parser = argparse.ArgumentParser(description="Step7 Scanpy processing with flexible inputs")
    parser.add_argument(
        '--inputMatrixDir', '-i', required=True,
        help="10X filtered matrix parent directory, e.g. 04-Matrix/filtered/"
    )
    parser.add_argument(
        '--inputfileBarcodelist', '-b', required=True,
        help="文件路径，包含barcode列表，每行一个barcode（1起始索引）"
    )
    parser.add_argument(
        '--splitCB', '-c', nargs='+', required=True,
        metavar='SPEC',
        help="索引规范列表，支持：\n"
             "1. 范围格式：1-5（包含起始和结束）\n"
             "2. 单个索引：3\n"
             "3. 等差数列：[1,89,8]（起始,结束,步长，1-based）"
    )
    parser.add_argument(
        '--splitSample', '-s', nargs='+', required=True,
        metavar='SAMPLE',
        help="与splitCB数量对应的样本名称列表，顺序需严格匹配"
    )
    return parser.parse_args()

def parse_index_spec(spec, total_length):
    """
    解析索引规范，支持三种格式：
    1. 范围格式：1-5 → 转换为1-based的[1,2,3,4,5]
    2. 单个索引：3 → 转换为[3]
    3. 等差数列：[1,89,8] → 1-based等差数列，包含起始和结束
    返回0-based索引列表（内部使用）
    """
    indices = []
    
    # 处理等差数列格式 [start,end,step]
    if spec.startswith('[') and spec.endswith(']'):
        try:
            content = spec[1:-1].strip()
            start, end, step = map(int, content.split(','))
            # 校验步长
            if step <= 0:
                raise ValueError(f"步长必须为正整数，当前值：{step}")
            # 处理顺序问题，确保start <= end
            if start > end:
                start, end = end, start
            # 生成1-based索引并转换为0-based
            current = start
            while current <= end:
                indices.append(current - 1)  # 转换为0-based
                current += step
        except Exception as e:
            raise ValueError(f"无效的等差数列格式 '{spec}': {e}")
        return indices
    
    # 处理范围或单个索引格式
    parts = spec.split(',')
    for part in parts:
        part = part.strip()
        if '-' in part:
            # 范围格式
            try:
                a, b = map(int, part.split('-'))
                if a < 1 or b < a or b > total_length:
                    raise ValueError
                # 生成1-based范围并转换为0-based
                indices.extend(range(a-1, b))  # range是左闭右开，b不转换因为原格式包含end
            except:
                raise ValueError(f"无效的范围格式 '{part}'，应为1-based的start-end")
        else:
            # 单个索引
            try:
                idx = int(part)
                if idx < 1 or idx > total_length:
                    raise ValueError
                indices.append(idx - 1)  # 转换为0-based
            except:
                raise ValueError(f"无效的单个索引 '{part}'，应为1-based整数")
    
    # 去重并保持顺序
    seen = set()
    return [x for x in indices if x not in seen and not seen.add(x)]

def main():
    args = parse_args()

    # 校验参数数量一致性
    if len(args.splitCB) != len(args.splitSample):
        raise ValueError(f"splitCB({len(args.splitCB)})与splitSample({len(args.splitSample)})数量不匹配")

    # 输入目录校验
    if not os.path.isdir(args.inputMatrixDir):
        raise FileNotFoundError(f"输入矩阵目录不存在: {args.inputMatrixDir}")
    
    # 加载barcode列表（1-based索引对应文件中的行号）
    with open(args.inputfileBarcodelist) as f:
        barcodes = [line.strip() for line in f if line.strip()]
    total_barcodes = len(barcodes)
    if total_barcodes == 0:
        raise ValueError(f"Barcode列表文件为空: {args.inputfileBarcodelist}")

    # 解析每个splitCB规范为0-based索引列表
    sample_indices = []
    for spec in args.splitCB:
        indices = parse_index_spec(spec, total_barcodes)
        # 校验索引范围
        for idx in indices:
            if idx < 0 or idx >= total_barcodes:
                raise ValueError(f"索引 {idx+1} 超出barcode列表范围（1-{total_barcodes}）")
        sample_indices.append(indices)

    # 构建barcode到sample的映射（使用1-based原始索引匹配barcodes列表）
    barcode_to_sample = {}
    for indices, sample in zip(sample_indices, args.splitSample):
        for idx in indices:
            bc = barcodes[idx]  # 通过0-based索引获取barcode
            if bc in barcode_to_sample:
                raise ValueError(f"Barcode {bc} 被分配到多个样本")
            barcode_to_sample[bc] = sample

    # 处理矩阵文件
    lib_dirs = glob(os.path.join(args.inputMatrixDir, "*_filtered_feature_bc_matrix"))
    if not lib_dirs:
        raise FileNotFoundError(f"未找到过滤矩阵目录，路径模式：{args.inputMatrixDir}/*_filtered_feature_bc_matrix")
    
    metadata_list = []
    adata_list = []
    for lib_dir in lib_dirs:
        lib_name = os.path.basename(lib_dir).replace("_filtered_feature_bc_matrix", "")
        print(f"处理文库: {lib_name}")
        adata = sc.read_10x_mtx(lib_dir, var_names='gene_symbols', cache=True)
        adata.obs['Library'] = lib_name
        metadata_list.append({**compute_library_stats(adata), 'Library': lib_name})
        adata_list.append(adata)
    
    # 保存单文库元数据
    pd.DataFrame(metadata_list).to_excel("05-Combined/3-PerLib_Metadata.xlsx", index=False)
    print("单文库元数据已保存")

    # 合并所有文库
    merged_adata = sc.concat(
        adata_list, join='outer', label='Library',
        keys=[d.obs['Library'].iloc[0] for d in adata_list], index_unique=None
    )
    sc.pp.calculate_qc_metrics(merged_adata, inplace=True)
    pd.DataFrame([compute_library_stats(merged_adata)]).to_excel("05-Combined/3-Merge_Metadata.xlsx", index=False)
    merged_adata.write("05-Combined/Merge.h5ad")
    print("合并文库及元数据已保存")

    # 分配样本标签
    merged_adata.obs['Barcode2'] = merged_adata.obs_names.str.split('_').str[1]  # 提取第二段barcode
    merged_adata.obs['Sample'] = merged_adata.obs['Barcode2'].map(barcode_to_sample)
    
    # 检查未匹配的barcode
    if merged_adata.obs['Sample'].isnull().any():
        unmatched = merged_adata.obs[merged_adata.obs['Sample'].isnull()]['Barcode2'].unique()
        print(f"警告：{len(unmatched)} 个Barcode未匹配到样本: {unmatched[:5]}...")

    # 按样本拆分并保存
    for sample in args.splitSample:
        sub_adata = merged_adata[merged_adata.obs['Sample'] == sample].copy()
        sub_adata.write(f"06-CB2_Sample/{sample}.h5ad")
        sub_metadata = compute_library_stats(sub_adata)
        sub_metadata['Sample'] = sample  # 添加样本标识
        pd.DataFrame([sub_metadata]).to_excel(f"06-CB2_Sample/4-{sample}_Metadata.xlsx", index=False)
        print(f"样本 {sample} 数据及元数据已保存")

    print("所有处理完成")

if __name__ == "__main__":
    main()