"""
SampleID_split.py

用法示例：
python SampleID_split.py \
  --inputfileCellReads 文件1 文件2 ... \
  --inputfileBarcodelist Barcode2.list \
  --inputfileCellSummaryFilter 2-Cell_Summary_withSampleID.xls \
  --splitCB [1,89,8] [2,90,8] [3,91,8] [4,92,8] [5,93,8] [6,94,8] [7,95,8] [8,96,8] \
  --splitSample A B C D E F G H \
  --inputfileCellMatrix MatrixA MatrixB MatrixC MatrixD MatrixE MatrixF MatrixG MatrixH \
  --outputfileCellSummaryRaw 2-Cell_Summary_withSampleID_AllCB.xls \
  --outputfileSampleSummary 5-Sample_Summary.xls

脚本原理（关键改进）：
1. 支持等差数列格式的splitCB参数（如[1,89,8]表示1-based索引，步长8的序列）
2. 改进索引解析逻辑，兼容范围/单个索引/等差数列三种格式
3. 优化样本分配逻辑，通过索引列表匹配Barcode所属样本
"""

import argparse
import pandas as pd
import numpy as np
import scanpy as sc
from pathlib import Path  # 新增：用于路径处理
import os

# ------------------ 函数定义 ------------------
def parse_index_spec(spec, total_length):
    """
    解析索引规范（兼容三种格式）：
    1. 等差数列：[start,end,step] 如[1,89,8]（1-based，包含起止，步长8）
    2. 范围格式：1-5（1-based，包含起止）
    3. 单个索引：3（1-based）
    
    返回：0-based索引列表（内部使用）
    """
    indices = []
    
    # 处理等差数列格式 [start,end,step]
    if spec.startswith('[') and spec.endswith(']'):
        try:
            content = spec[1:-1].strip()
            start, end, step = map(int, content.split(','))
            # 校验步长必须为正整数
            if step <= 0:
                raise ValueError(f"步长必须>0，当前值：{step}")
            # 确保start <= end
            if start > end:
                start, end = end, start
            # 生成1-based索引并转换为0-based
            current = start
            while current <= end:
                indices.append(current - 1)  # 转换为0-based索引
                current += step
        except Exception as e:
            raise ValueError(f"无效等差数列格式 '{spec}': {e}")
        return indices
    
    # 处理范围或单个索引格式
    parts = spec.split(',')
    for part in parts:
        part = part.strip()
        if '-' in part:
            # 范围格式：1-5 → 1-based
            try:
                a, b = map(int, part.split('-'))
                if a < 1 or b < a or b > total_length:
                    raise ValueError
                indices.extend(range(a-1, b))  # 转换为0-based范围（左闭右开，兼容原格式包含end）
            except:
                raise ValueError(f"无效范围格式 '{part}'，应为1-based的start-end")
        else:
            # 单个索引：3 → 1-based
            try:
                idx = int(part)
                if idx < 1 or idx > total_length:
                    raise ValueError
                indices.append(idx - 1)  # 转换为0-based索引
            except:
                raise ValueError(f"无效单个索引 '{part}'，应为1-based整数")
    
    # 去重并保持顺序
    seen = set()
    return [x for x in indices if x not in seen and not seen.add(x)]

# ------------------ 主流程 ------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='SampleID split and summary (支持等差数列索引)')
    parser.add_argument('--inputfileCellReads', nargs='+', help='CellReads.stats 文件列表', required=True)
    parser.add_argument('--inputfileBarcodelist', help='Barcode 库文件（每行一个barcode，1-based索引对应行号）', required=True)
    parser.add_argument('--inputfileCellSummaryFilter', help='过滤后带SampleID的CellSummary文件', required=True)
    parser.add_argument('--splitCB', nargs='+', required=True,
                        metavar='SPEC',
                        help='索引规范列表，支持：\n'
                             '1. 等差数列：[1,89,8]（起始,结束,步长，1-based）\n'
                             '2. 范围格式：1-5\n'
                             '3. 单个索引：3')
    parser.add_argument('--splitSample', nargs='+', required=True,
                        metavar='SAMPLE',
                        help='与splitCB数量对应的样本名称列表（顺序严格匹配）')
    parser.add_argument('--inputfileCellMatrix', nargs='+', help='filtered_feature_bc_matrix目录列表', required=True)
    parser.add_argument('--outputfileCellSummaryRaw', help='输出带SampleID的原始表', required=True)
    parser.add_argument('--outputfileSampleSummary', help='输出样本汇总表', required=True)
    args = parser.parse_args()

    # ====================== 新增路径处理逻辑 ======================
    # 定义矩阵文件的基路径（所有矩阵目录均位于此子目录下）
    MATRIX_BASE_DIR = Path("07-CB2_Matrix")

    # 校验基路径是否存在（可选，根据实际需求决定是否强制创建）
    if not MATRIX_BASE_DIR.exists():
        print(f"警告: 基路径 {MATRIX_BASE_DIR} 不存在，脚本将尝试自动创建")
        MATRIX_BASE_DIR.mkdir(parents=True, exist_ok=True)

    # 自动为每个矩阵简称添加基路径前缀
    args.inputfileCellMatrix = [
        str(MATRIX_BASE_DIR / matrix_name)  # 使用Path拼接路径，自动处理分隔符
        for matrix_name in args.inputfileCellMatrix
    ]
    # ==============================================================

    # 校验参数数量一致性
    if len(args.splitCB) != len(args.splitSample):
        raise ValueError(f"splitCB({len(args.splitCB)})与splitSample({len(args.splitSample)})数量不匹配")

    # 读取Barcode列表（1-based索引对应文件中的行号）
    with open(args.inputfileBarcodelist) as f:
        barcode_list = [line.strip() for line in f if line.strip()]
    total_length = len(barcode_list)
    if total_length == 0:
        raise ValueError(f"Barcode列表文件为空: {args.inputfileBarcodelist}")

    # 解析splitCB规范为0-based索引列表
    sample_indices = []
    for spec in args.splitCB:
        indices = parse_index_spec(spec, total_length)
        # 校验索引范围
        for idx in indices:
            if idx < 0 or idx >= total_length:
                raise ValueError(f"索引 {idx+1} 超出Barcode列表范围（1-{total_length}）")
        sample_indices.append(indices)

    # 建立barcode到样本的映射（通过0-based索引获取barcode）
    barcode_to_sample = {}
    for indices, sample in zip(sample_indices, args.splitSample):
        for idx in indices:
            bc = barcode_list[idx]  # 通过0-based索引获取barcode
            if bc in barcode_to_sample:
                raise ValueError(f"Barcode {bc} 被分配到多个样本")
            barcode_to_sample[bc] = sample

    # 步骤1+2：处理并合并CellReads.stats（原逻辑不变）
    df_list = []
    for fn in args.inputfileCellReads:
        df = pd.read_csv(fn, sep='\t', header=0, skiprows=[1], dtype={'CB': str})
        df_list.append(df)
    df_all = pd.concat(df_list, axis=0, ignore_index=True)
    df_all.sort_values(by=df_all.columns[16], ascending=False, inplace=True)

    # 步骤3：提取第二段barcode并分配SampleID（关键修改：通过索引列表匹配）
    def assign_sample(cb):
        """
        从CB中提取第二段barcode，根据索引列表分配样本
        """
        segs = cb.split('_')
        if len(segs) < 3:
            return 'Unknown'
        mid_bc = segs[1]
        # 检查barcode是否在映射表中
        return barcode_to_sample.get(mid_bc, 'Unknown')

    df_all['SampleID'] = df_all['CB'].apply(assign_sample)
    df_all.to_csv(args.outputfileCellSummaryRaw, sep='\t', index=False)
    print(f"已生成带SampleID的原始表: {args.outputfileCellSummaryRaw}")

    # 步骤4-6：计算各项指标（原逻辑不变，仅修改注释）
    # [中间计算代码保持不变，仅调整注释说明参数格式]

    # 步骤7：整合指标（原逻辑不变）
    # [后续代码与原脚本一致，省略重复部分...]

    # （以下为原脚本的指标计算和输出部分，保持不变）
    # 步骤4：基于CellReads原始表计算整体mapping指标
    grp = df_all.groupby('SampleID')
    reads_valid = grp['cbMatch'].sum()
    uniq_map = grp['genomeU'].sum()
    uniq_map_frac = (uniq_map / reads_valid * 100).round(2).astype(str)
    multi_map = grp['genomeM'].sum()
    multi_map_frac = (multi_map / reads_valid * 100).round(2).astype(str)
    exonic = grp['exonic'].sum()
    intronic = grp['intronic'].sum()
    antisense = (grp['exonicAS'].sum() + grp['intronicAS'].sum())
    sum6_7 = uniq_map + multi_map
    sum10_13 = exonic + intronic + grp['exonicAS'].sum() + grp['intronicAS'].sum()
    intergenic = sum6_7 - sum10_13
    assigned_sum = sum6_7
    exonic_frac = (exonic / assigned_sum * 100).round(2).astype(str)
    intronic_frac = (intronic / assigned_sum * 100).round(2).astype(str)
    intergenic_frac = (intergenic / assigned_sum * 100).round(2).astype(str)
    antisense_frac = (antisense / assigned_sum * 100).round(2).astype(str)

    # 步骤5：处理过滤后CellSummary表
    df2 = pd.read_csv(args.inputfileCellSummaryFilter, sep='\t', header=0, dtype={'CB': str})
    g2 = df2.groupby('SampleID')
    est_cells = g2.size().astype(int)
    uniq_reads_cells = g2['countedU'].sum().astype(int)
    mean_genic = g2['countedU'].mean().round(0).astype(int)
    mean_umi = g2['nUMIunique'].mean().round(0).astype(int)
    mean_gene = g2['nGenesUnique'].mean().round(0).astype(int)
    med_genic = g2['countedU'].median().round(0).astype(int)
    med_umi = g2['nUMIunique'].median().round(0).astype(int)
    med_gene = g2['nGenesUnique'].median().round(0).astype(int)
    countedU_sum_filt = g2['countedU'].sum()
    nUMIunique_sum_filt = g2['nUMIunique'].sum()
    seq_sat = (1 - nUMIunique_sum_filt / countedU_sum_filt) * 100
    seq_sat = seq_sat.round(2).astype(str)

    # # 步骤6：计算基因覆盖数
    # genes_covered = {}
    # for matrix_dir, sample in zip(args.inputfileCellMatrix, args.splitSample):
    #     adata = sc.read_10x_mtx(matrix_dir, cache=False)
    #     X = adata.X
    #     try:
    #         gene_nonzero = np.array((X.sum(axis=0) > 0)).reshape(-1)
    #     except:
    #         gene_nonzero = (X.sum(axis=0) > 0)
    #     genes_covered[sample] = int(gene_nonzero.sum())
    # 步骤6：计算基因覆盖数（修改点：使用处理后的路径）
    genes_covered = {}
    for matrix_dir, sample in zip(args.inputfileCellMatrix, args.splitSample):
        # matrix_dir 已自动添加 "07-CB2_Matrix/" 前缀，例如："07-CB2_Matrix/A_filtered_feature_bc_matrix"
        try:
            adata = sc.read_10x_mtx(matrix_dir, cache=False)  # 直接使用完整路径读取
            X = adata.X
            gene_nonzero = np.array((X.sum(axis=0) > 0)).reshape(-1)
            genes_covered[sample] = int(gene_nonzero.sum())
        except Exception as e:
            raise ValueError(f"读取矩阵文件 {matrix_dir} 失败: {e}")

    # 处理Mean Reads per Cell（避免NaN/inf）
    lib_mean = pd.Series(dtype=int)
    for sample in args.splitSample:
        rv = reads_valid.get(sample, 0)
        ec = est_cells.get(sample, 1)  # 避免除零
        lib_mean[sample] = round(rv / ec) if ec != 0 else 0

    # 整合指标（保持原顺序）
    metrics = [
        ('Reads With Valid Barcodes Num', reads_valid),
        ('Sequencing Saturation', seq_sat),
        ('Uniquely Mapped Reads', uniq_map),
        ('Uniquely Mapped Reads Fraction', uniq_map_frac),
        ('Multi-Mapped Reads', multi_map),
        ('Multi-Mapped Reads Fraction', multi_map_frac),
        ('Mapped Reads Assigned To Exonic Regions', exonic),
        ('Mapped Reads Assigned To Intronic Regions', intronic),
        ('Mapped Reads Assigned Antisense To Gene', antisense),
        ('Mapped Reads Assigned To Intergenic Regions', intergenic),
        ('Mapped Reads Assigned Sum', assigned_sum),
        ('Mapped Reads Assigned To Exonic Regions Fraction', exonic_frac),
        ('Mapped Reads Assigned To Intronic Regions Fraction', intronic_frac),
        ('Mapped Reads Assigned To Intergenic Regions Fraction', intergenic_frac),
        ('Mapped Reads Assigned Antisense To Gene Fraction', antisense_frac),
        ('Estimated Number of Cells', est_cells),
        ('Unique Reads in Cells Mapped to Gene', uniq_reads_cells),
        ('Mean Reads per Cell (based on Reads With Valid Barcodes Num)', lib_mean),
        ('Mean GenicReads per Cell', mean_genic),
        ('Mean UMI per Cell', mean_umi),
        ('Mean Gene per Cell', mean_gene),
        ('Median GenicReads per Cell', med_genic),
        ('Median UMI per Cell', med_umi),
        ('Median Gene per Cell', med_gene),
        ('Total Gene Detected', pd.Series(genes_covered))
    ]

    # 构建并保存结果
    summary_df = pd.DataFrame({sam: {name: vals.get(sam, np.nan) for name, vals in metrics} for sam in args.splitSample})
    # summary_df = summary_df.T
    summary_df.to_csv(args.outputfileSampleSummary, sep='\t', index=True, header=True)
    print(f"已生成样本汇总表: {args.outputfileSampleSummary}")