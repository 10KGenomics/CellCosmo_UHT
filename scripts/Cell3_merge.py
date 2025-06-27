"""
Cell3_merge.py: 合并和汇总单细胞测序子文库分析结果脚本（含第28步Sequencing Saturation计算），并新增Total Genes Covered指标

功能:
  1. 从 Results_Summary.xls 读取所有子文库的统计数据，计算汇总指标（Reads With Valid Barcodes, Mapping, 区域分配等）。
  2. 根据参数 SplitSample 指定的批次，分别计算 Number of Reads 的均值并汇总；如果 SplitSample=False，则统一批次，求和后取平均。
  3. 计算新的 Reads With Valid Barcodes（百分比），保留两位小数。
  4. 计算 Q30 Bases in CB+UMI / Q30 Bases in RNA read 的均值，保留两位小数。
  5. 计算 Uniquely / Multi-Mapped Reads 总和及其百分比（保留两位小数）。
  6. 计算 Mapped Reads Assigned 到 Exonic/Intronic/Intergenic/Antisense 的总和及其分数（保留两位小数）。
  7. 统计 Estimated Number of Cells、Unique Reads in Cells Mapped to Gene 总和，Fraction of Unique Reads 平均值（保留两位小数）。
  8. 读取并处理每个子文库的 CellReads.stats 文件，生成 Cell_Summary.xls。
  9. 计算 Cell_Summary.xls 第15/17/18列的均值和中位数。
 10. 计算 Sequencing Saturation 并插入至汇总指标中。
 11. 新增 Mean Reads per Cell 指标并插入至汇总指标中。
 12. 【新增】读取所有10X矩阵目录，合并矩阵，统计基因表达>0的基因数，作为 Total Genes Covered 指标。
 13. 最终将所有指标写入 Cell3_Merge_Summary.xls。

使用示例:
  python Cell3_merge.py \
    --SplitSample 16to1 8to1 \
    --inputfileSummary 1-Results_Summary.xls \
    --inputfileCellReads $(for i in $(ls 02-STARsolo/*-Solo.out/GeneFull_Ex50pAS/CellReads.stats);do echo " $i";done) \
    --inputfileMatrixDir $(for i in $(ls -d 04-Matrix/filtered/*_filtered_feature_bc_matrix);do echo " $i";done) \
    --outputfile1 2-Cell_Summary.xls \
    --outputfile2 2-Cell3_Merge_Summary.xls
"""

import os
import glob
import argparse
import pandas as pd
import numpy as np
import scipy.io
from scipy import sparse


def parse_args():
    parser = argparse.ArgumentParser(description='合并并汇总单细胞测序子文库结果，新增Total Genes Covered指标')
    parser.add_argument('--SplitSample', nargs='+', required=True,
                        help='批次名称列表，如：16to1 8to1；或 False')
    parser.add_argument('--inputfileSummary', required=True,
                        help='Results_Summary.xls 文件路径')
    parser.add_argument('--inputfileCellReads', nargs='+', required=True,
                        help='各子文库 CellReads.stats 文件列表')
    parser.add_argument('--inputfileMatrixDir', nargs='+', required=True,
                        help='10X过滤后矩阵目录列表，每个目录应包含 matrix.mtx(.gz), features.tsv(.gz), barcodes.tsv(.gz)')
    parser.add_argument('--outputfile1', required=True,
                        help='输出合并后细胞文件 Cell_Summary.xls')
    parser.add_argument('--outputfile2', required=True,
                        help='输出最终汇总指标 Cell3_Merge_Summary.xls')
    return parser.parse_args()


def summarize_library_metrics(df, splits):
    """
    基于 Results_Summary DataFrame 和 SplitSample 参数，计算库级汇总指标。
    返回一个按插入顺序保存指标的 dict（不含Sequencing Saturation和Total Genes Covered）。
    """
    metrics = {}

    # 1. 计算新的 Number of Reads
    num_reads = df.loc['Number of Reads']
    if len(splits) == 1 and splits[0].lower() == 'false':
        new_nr = num_reads.sum() / len(num_reads)
    else:
        vals = []
        for batch in splits:
            cols = [c for c in df.columns if c.startswith(batch + '_')]
            if not cols:
                raise ValueError(f"未找到以'{batch}_'开头的列")
            vals.append(num_reads[cols].mean())
        new_nr = sum(vals)
    metrics['Number of Reads'] = float(new_nr)

    # 2. Reads With Valid Barcodes Num 总和
    sum_rwvbn = int(df.loc['Reads With Valid Barcodes Num'].sum())
    metrics['Reads With Valid Barcodes Num'] = sum_rwvbn

    # 3. 新 Reads With Valid Barcodes 百分比
    rwvb_pct = round(sum_rwvbn / new_nr * 100, 2)
    metrics['Reads With Valid Barcodes'] = rwvb_pct

    # 4. Q30 均值
    metrics['Q30 Bases in CB+UMI'] = round(df.loc['Q30 Bases in CB+UMI'].mean(), 2)
    metrics['Q30 Bases in RNA read'] = round(df.loc['Q30 Bases in RNA read'].mean(), 2)

    # 5-6. Mapped Reads 总和
    umap = int(df.loc['Uniquely Mapped Reads'].sum())
    mmap = int(df.loc['Multi-Mapped Reads'].sum())
    metrics['Uniquely Mapped Reads'] = umap
    metrics['Multi-Mapped Reads'] = mmap

    # 7-8. Mapped Reads Fraction
    metrics['Uniquely Mapped Reads Fraction'] = round(umap / sum_rwvbn * 100, 2)
    metrics['Multi-Mapped Reads Fraction'] = round(mmap / sum_rwvbn * 100, 2)

    # 9-12. 区域分配总和
    exonic = int(df.loc['Mapped Reads Assigned To Exonic Regions'].sum())
    intronic = int(df.loc['Mapped Reads Assigned To Intronic Regions'].sum())
    intergenic = int(df.loc['Mapped Reads Assigned To Intergenic Regions'].sum())
    antisense = int(df.loc['Mapped Reads Assigned Antisense To Gene'].sum())
    metrics['Mapped Reads Assigned To Exonic Regions'] = exonic
    metrics['Mapped Reads Assigned To Intronic Regions'] = intronic
    metrics['Mapped Reads Assigned To Intergenic Regions'] = intergenic
    metrics['Mapped Reads Assigned Antisense To Gene'] = antisense

    # 13. Mapped Reads Assigned Sum
    total_map = exonic + intronic + intergenic + antisense
    metrics['Mapped Reads Assigned Sum'] = total_map

    # 14-17. 区域分配 Fraction
    metrics['Mapped Reads Assigned To Exonic Regions Fraction'] = round(exonic / total_map * 100, 2)
    metrics['Mapped Reads Assigned To Intronic Regions Fraction'] = round(intronic / total_map * 100, 2)
    metrics['Mapped Reads Assigned To Intergenic Regions Fraction'] = round(intergenic / total_map * 100, 2)
    metrics['Mapped Reads Assigned Antisense To Gene Fraction'] = round(antisense / total_map * 100, 2)

    # 18-20. 细胞相关库级指标
    metrics['Estimated Number of Cells'] = int(df.loc['Estimated Number of Cells'].sum())
    metrics['Unique Reads in Cells Mapped to Gene'] = int(
        df.loc['Unique Reads in Cells Mapped to Gene'].sum()
    )
    metrics['Fraction of Unique Reads in Cells'] = round(
        df.loc['Fraction of Unique Reads in Cells'].mean(), 2
    )

    return metrics


def process_stats(filepath, est_cells):
    """
    处理单库 CellReads.stats: 去除前两行 -> 排序第17列 -> 保留前 est_cells 行
    返回生成的 tmp3 文件路径
    """
    tmp1 = filepath + '.tmp1'
    tmp2 = filepath + '.tmp2'
    tmp3 = filepath + '.tmp3'

    # 去除前两行无效信息
    with open(filepath) as fin, open(tmp1, 'w') as fout:
        fout.writelines(fin.readlines()[2:])

    # 排序并写 tmp2
    df = pd.read_csv(tmp1, sep='\t', header=None)
    df_sorted = df.sort_values(by=16, ascending=False)
    df_sorted.to_csv(tmp2, sep='\t', header=False, index=False)

    # 截取前 est_cells 行，写 tmp3
    df_top = df_sorted.head(est_cells)
    df_top.to_csv(tmp3, sep='\t', header=False, index=False)
    return tmp3


def compute_total_genes_covered(matrix_dirs):
    """
    读取所有10X filtered_feature_bc_matrix目录，合并矩阵后统计基因表达>0的基因数
    返回 Total Genes Covered (整数)
    """
    mats = []
    gene_index = None
    for md in matrix_dirs:
        # 查找 matrix 文件（支持 .mtx 或 .mtx.gz）
        mtx_candidates = glob.glob(os.path.join(md, 'matrix.mtx*'))
        if not mtx_candidates:
            raise ValueError(f"在目录 {md} 中未找到 matrix.mtx 文件")
        mtx_path = mtx_candidates[0]
        # 读取稀疏矩阵（Matrix Market 格式，支持 gzip）
        mat = scipy.io.mmread(mtx_path).tocsc()

        # 查找 features 文件（支持 features.tsv 或 genes.tsv，及 .gz 压缩）
        feat_candidates = glob.glob(os.path.join(md, 'features.tsv*')) + glob.glob(os.path.join(md, 'genes.tsv*'))
        if not feat_candidates:
            raise ValueError(f"在目录 {md} 中未找到 features.tsv 或 genes.tsv 文件")
        feats_path = feat_candidates[0]
        # 自动根据后缀处理 gzip
        compression = 'gzip' if feats_path.endswith('.gz') else None
        feats = pd.read_csv(feats_path, sep='\t', header=None, compression=compression)
        genes = feats[1].values

        # 初始化或校验基因索引一致性
        if gene_index is None:
            gene_index = genes
        else:
            if not np.array_equal(gene_index, genes):
                raise ValueError(f"基因列表在 {md} 中与首个矩阵不一致")

        mats.append(mat)

    # 合并所有细胞矩阵（按列拼接）
    combined = sparse.hstack(mats)
    # 计算每个基因在所有细胞中的表达总和
    gene_sums = np.array(combined.sum(axis=1)).flatten()
    # 统计表达量>0的基因数
    total_genes = int((gene_sums > 0).sum())
    return total_genes


def main():
    args = parse_args()
    # 读取 Results_Summary.xls
    summary_df = pd.read_csv(args.inputfileSummary, sep='\t', index_col=0)

    # 计算库级指标（不含Sequencing Saturation, Mean Reads per Cell, Total Genes Covered）
    metrics = summarize_library_metrics(summary_df, args.SplitSample)

    # 处理并合并所有子库的 tmp3 文件，生成 Cell_Summary.xls
    tmp3_paths = []
    for stat in args.inputfileCellReads:
        lib = os.path.basename(os.path.dirname(os.path.dirname(stat))).replace('-Solo.out','')
        est = int(summary_df.loc['Estimated Number of Cells', lib])
        tmp3_paths.append(process_stats(stat, est))

    combined = pd.concat(
        [pd.read_csv(p, sep='\t', header=None) for p in tmp3_paths],
        ignore_index=True
    )
    combined.to_csv(args.outputfile1, sep='\t', header=False, index=False)

    # # 计算细胞级均值和中位数
    # mean_r = int(combined[14].mean())
    # mean_u = int(combined[16].mean())
    # mean_g = int(combined[17].mean())
    # med_r  = int(combined[14].median())
    # med_u  = int(combined[16].median())
    # med_g  = int(combined[17].median())
    # ——更正为：先四舍五入，再转换为 int——  
    mean_r = int(round(combined[14].mean())) 
    mean_u = int(round(combined[16].mean()))
    mean_g = int(round(combined[17].mean()))
    med_r  = int(round(combined[14].median()))
    med_u  = int(round(combined[16].median()))
    med_g  = int(round(combined[17].median()))

    # 计算 Sequencing Saturation
    sum_r = combined[14].sum()
    sum_u = combined[16].sum()
    sat   = round((1 - sum_u / sum_r) * 100, 2)

    # 在库级指标中插入 Sequencing Saturation
    items = list(metrics.items())
    idx = next(i for i, (k, _) in enumerate(items) if k == 'Reads With Valid Barcodes')
    items.insert(idx + 1, ('Sequencing Saturation', sat))

    # 新增 Mean Reads per Cell
    lib_mean = int(metrics['Reads With Valid Barcodes Num'] / metrics['Estimated Number of Cells'])
    i2 = next(i for i,(k,_) in enumerate(items) if k=='Fraction of Unique Reads in Cells')
    items.insert(i2+1, ('Mean Reads per Cell(based on Reads With Valid Barcodes Num)', lib_mean))

    # 新增 Mean Reads per Cell
    lib_mean_raw = int(metrics['Number of Reads'] / metrics['Estimated Number of Cells'])
    i2 = next(i for i,(k,_) in enumerate(items) if k=='Fraction of Unique Reads in Cells')
    items.insert(i2+1, ('Mean Reads per Cell(based on Number of RawReads)', lib_mean_raw))    

    # 新增 Total Gene Detected
    total_genes = compute_total_genes_covered(args.inputfileMatrixDir)
    items.append(('Total Gene Detected', total_genes))

    # 组装并写出最终汇总表
    rows = [{'Sample': k, 'Merged': v} for k, v in items]
    rows.extend([
        {'Sample': 'Mean GenicReads per Cell', 'Merged': mean_r},
        {'Sample': 'Mean UMI per Cell', 'Merged': mean_u},
        {'Sample': 'Mean Gene per Cell', 'Merged': mean_g},
        {'Sample': 'Median GenicReads per Cell', 'Merged': med_r},
        {'Sample': 'Median UMI per Cell', 'Merged': med_u},
        {'Sample': 'Median Gene per Cell', 'Merged': med_g},
    ])
    final_df = pd.DataFrame(rows)
    final_df.to_csv(args.outputfile2, sep='\t', index=False)
    #final_df.to_csv(args.outputfile2, sep='\t', index=False, float_format='%.2f')  # 保留两位小数，整数自动处理

    print(f"合并细胞数据写入: {args.outputfile1}")
    print(f"最终汇总指标写入: {args.outputfile2}")

if __name__ == '__main__':
    main()