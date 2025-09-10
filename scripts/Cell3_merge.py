"""
Cell3_merge.py: 合并和汇总单细胞测序子文库分析结果脚本（含Sequencing Saturation计算）
新增功能：从splitcode生成的summary文件（my_summary.txt）提取总测序reads数，新增Total Gene Detected指标

功能概述:
  1. 从splitcode的summary文件（JSON格式）提取总测序读数（n_reads_max）作为Number of Reads
  2. 从Results_Summary.xls读取子文库统计数据，计算汇总指标（有效Barcode、比对率等）
  3. 处理CellReads.stats文件，生成细胞级汇总表（Cell_Summary.xls）
  4. 计算测序饱和度（Sequencing Saturation）、细胞平均reads数等关键指标
  5. 合并10X矩阵，统计检测到的总基因数（Total Gene Detected）
  6. 输出最终汇总指标至Cell3_Merge_Summary.xls

使用示例:
  python Cell3_merge.py \
    --SplitSample pbmc \
    --inputfileSummary 05-Combined/1-Results_Summary.xls \
    --inputfileSplitSummary my_summary.txt \
    --inputfileCellReads $(ls 02-STARsolo/star_outs_lib*/Solo.out/GeneFull_Ex50pAS/CellReads.stats) \
    --inputfileMatrixDir $(ls -d 04-CB3_Matrix/filtered/*_filtered_feature_bc_matrix) \
    --outputfile1 05-Combined/2-Cell_Summary.xls \
    --outputfile2 05-Combined/2-Cell3_Merge_Summary.xls
"""

import os
import re
import glob
import json
import argparse
import pandas as pd
import numpy as np
import scipy.io
from scipy import sparse
from pathlib import Path  # 用于路径处理，更直观


def parse_args():
    """解析命令行参数，定义输入输出文件路径及分析参数"""
    parser = argparse.ArgumentParser(description='合并单细胞子文库结果并计算汇总指标（含splitcode summary解析）')
    # 批次名称参数（如样本名或False表示不拆分）
    parser.add_argument('--SplitSample', nargs='+', required=True,
                        help='批次名称列表（如pbmc）或False（不拆分批次），例：pbmc 或 False')
    # 原有子文库统计文件（Results_Summary.xls）
    parser.add_argument('--inputfileSummary', required=True,
                        help='子文库统计结果文件（如1-Results_Summary.xls），包含各子库基础指标')
    # 新增：splitcode生成的summary文件（my_summary.txt，JSON格式），用于提取总reads数
    parser.add_argument('--inputfileSplitSummary', required=True,
                        help='splitcode输出的JSON格式汇总文件（如my_summary.txt），需包含n_reads_max字段')
    # 各子库的CellReads.stats文件列表
    parser.add_argument('--inputfileCellReads', nargs='+', required=True,
                        help='所有子文库的CellReads.stats文件路径，用于细胞级指标计算')
    # 10X过滤后矩阵目录列表
    parser.add_argument('--inputfileMatrixDir', nargs='+', required=True,
                        help='10X过滤矩阵目录（含matrix.mtx、features.tsv、barcodes.tsv），用于统计总基因数')
    # 输出文件1：细胞级汇总表
    parser.add_argument('--outputfile1', required=True,
                        help='输出细胞级汇总结果（如2-Cell_Summary.xls）')
    # 输出文件2：最终汇总指标表
    parser.add_argument('--outputfile2', required=True,
                        help='输出最终合并的所有汇总指标（如2-Cell3_Merge_Summary.xls）')
    return parser.parse_args()


def extract_n_reads_max(split_summary_path):
    """
    从splitcode的JSON格式summary文件中提取n_reads_max的值
    
    参数:
        split_summary_path: splitcode输出的summary文件路径（如my_summary.txt）
    
    返回:
        int: n_reads_max的值（总测序reads数）
    
    异常:
        FileNotFoundError: 若文件不存在
        json.JSONDecodeError: 若文件不是有效的JSON格式
        ValueError: 若文件中不包含n_reads_max字段或字段值非整数
    """
    # 检查文件是否存在
    if not os.path.exists(split_summary_path):
        raise FileNotFoundError(f"splitcode summary文件不存在: {split_summary_path}")
    
    # 读取并解析JSON文件
    try:
        with open(split_summary_path, 'r') as f:
            summary_data = json.load(f)
    except json.JSONDecodeError as e:
        raise ValueError(f"JSON格式解析错误（{split_summary_path}）: {str(e)}")
    
    # 提取n_reads_max并验证
    n_reads_max = summary_data.get('n_reads_max')
    if n_reads_max is None:
        raise ValueError(f"文件{split_summary_path}中未找到'n_reads_max'字段")
    if not isinstance(n_reads_max, int):
        raise TypeError(f"n_reads_max必须是整数，实际为: {type(n_reads_max)}")
    
    return n_reads_max


def summarize_library_metrics(df, splits, n_reads_max):
    """
    基于Results_Summary数据和总reads数，计算库级汇总指标
    
    参数:
        df: Results_Summary.xls读取的DataFrame（行：指标，列：子库）
        splits: 批次名称列表（来自--SplitSample）
        n_reads_max: 从split_summary提取的总reads数（Number of Reads）
    
    返回:
        dict: 按顺序存储的库级汇总指标（不含测序饱和度、总基因数等后计算指标）
    """
    metrics = {}
    
    # 1. 总测序reads数（从split_summary提取）
    metrics['Number of Reads'] = n_reads_max
    
    # 2. 有效Barcode reads总数（各子库求和）
    sum_valid_barcode_reads = int(df.loc['Reads With Valid Barcodes Num'].sum())
    metrics['Reads With Valid Barcodes Num'] = sum_valid_barcode_reads
    
    # 3. 有效Barcode reads百分比（保留2位小数）
    valid_barcode_pct = round(sum_valid_barcode_reads / n_reads_max * 100, 2)
    metrics['Reads With Valid Barcodes'] = valid_barcode_pct
    
    # 4. Q30指标均值（CB+UMI和RNA read的Q30均值）
    # metrics['Q30 Bases in CB+UMI'] = round(df.loc['Q30 Bases in CB+UMI'].mean(), 2)
    metrics['Q30 Bases in RNA read'] = round(df.loc['Q30 Bases in RNA read'].mean(), 2)
    
    # 5. 比对reads总数（唯一比对+多重比对）
    unique_mapped = int(df.loc['Uniquely Mapped Reads'].sum())
    multi_mapped = int(df.loc['Multi-Mapped Reads'].sum())
    metrics['Uniquely Mapped Reads'] = unique_mapped
    metrics['Multi-Mapped Reads'] = multi_mapped
    
    # 6. 比对百分比（相对于有效Barcode reads）
    metrics['Uniquely Mapped Reads Fraction'] = round(unique_mapped / sum_valid_barcode_reads * 100, 2)
    metrics['Multi-Mapped Reads Fraction'] = round(multi_mapped / sum_valid_barcode_reads * 100, 2)
    
    # 7. 基因组区域分配总数（外显子、内含子、基因间、反义链）
    exonic = int(df.loc['Mapped Reads Assigned To Exonic Regions'].sum())
    intronic = int(df.loc['Mapped Reads Assigned To Intronic Regions'].sum())
    intergenic = int(df.loc['Mapped Reads Assigned To Intergenic Regions'].sum())
    antisense = int(df.loc['Mapped Reads Assigned Antisense To Gene'].sum())
    metrics['Mapped Reads Assigned To Exonic Regions'] = exonic
    metrics['Mapped Reads Assigned To Intronic Regions'] = intronic
    metrics['Mapped Reads Assigned To Intergenic Regions'] = intergenic
    metrics['Mapped Reads Assigned Antisense To Gene'] = antisense
    
    # 8. 区域分配总和及各区域占比
    total_region_mapped = exonic + intronic + intergenic + antisense
    metrics['Mapped Reads Assigned Sum'] = total_region_mapped
    metrics['Mapped Reads Assigned To Exonic Regions Fraction'] = round(exonic / total_region_mapped * 100, 2)
    metrics['Mapped Reads Assigned To Intronic Regions Fraction'] = round(intronic / total_region_mapped * 100, 2)
    metrics['Mapped Reads Assigned To Intergenic Regions Fraction'] = round(intergenic / total_region_mapped * 100, 2)
    metrics['Mapped Reads Assigned Antisense To Gene Fraction'] = round(antisense / total_region_mapped * 100, 2)
    
    # 9. 细胞相关指标（估计细胞数、细胞中唯一比对reads等）
    metrics['Estimated Number of Cells'] = int(df.loc['Estimated Number of Cells'].sum())
    metrics['Unique Reads in Cells Mapped to Gene'] = int(df.loc['Unique Reads in Cells Mapped to Gene'].sum())
    metrics['Fraction of Unique Reads in Cells'] = round(df.loc['Fraction of Unique Reads in Cells'].mean(), 2)
    
    return metrics


def process_cell_reads(stat_path, est_cells):
    """
    处理单个子库的CellReads.stats文件：过滤无效行、排序、截取前N个细胞（N=估计细胞数）
    
    参数:
        stat_path: CellReads.stats文件路径
        est_cells: 该子库的估计细胞数（从Results_Summary提取）
    
    返回:
        str: 处理后的临时文件路径（包含前N个细胞数据）
    
    步骤:
        1. 去除前2行无效信息
        2. 按第17列（UMI数）降序排序
        3. 保留前est_cells行（估计的细胞数）
        4. 清理中间临时文件，仅返回最终结果文件路径
    """
    # 定义临时文件路径（基于原文件路径，添加后缀区分）
    tmp_remove_header = f"{stat_path}.tmp_remove_header"  # 去除头两行的临时文件
    tmp_sorted = f"{stat_path}.tmp_sorted"  # 排序后的临时文件
    tmp_top_cells = f"{stat_path}.tmp_top_cells"  # 截取前N细胞的临时文件
    
    try:
        # 步骤1：去除前两行无效信息
        with open(stat_path, 'r') as fin, open(tmp_remove_header, 'w') as fout:
            fout.writelines(fin.readlines()[2:])  # 跳过前2行
        
        # 步骤2：按第17列（索引16，UMI数）降序排序
        df = pd.read_csv(tmp_remove_header, sep='\t', header=None)
        df_sorted = df.sort_values(by=16, ascending=False)  # 第17列（UMI）越大越可能是真实细胞
        df_sorted.to_csv(tmp_sorted, sep='\t', header=False, index=False)
        
        # 步骤3：截取前est_cells行（估计的细胞数）
        df_top = df_sorted.head(est_cells)
        df_top.to_csv(tmp_top_cells, sep='\t', header=False, index=False)
        
    finally:
        # 清理中间临时文件（仅保留最终结果文件）
        for tmp in [tmp_remove_header, tmp_sorted]:
            if os.path.exists(tmp):
                os.remove(tmp)
    
    return tmp_top_cells


def compute_total_genes_covered(matrix_dirs):
    """
    合并所有10X过滤矩阵，统计表达量>0的基因总数（Total Gene Detected）
    
    参数:
        matrix_dirs: 10X过滤矩阵目录列表（每个目录含matrix.mtx、features.tsv等）
    
    返回:
        int: 表达量>0的基因总数
    
    验证:
        确保所有矩阵的基因列表一致（避免批次间基因注释差异导致统计错误）
    """
    matrices = []  # 存储各子库的稀疏矩阵
    gene_list = None  # 用于验证所有子库基因列表一致性
    
    for dir_path in matrix_dirs:
        dir_path = Path(dir_path)  # 使用pathlib简化路径操作
        
        # 1. 查找矩阵文件（支持.mtx或.mtx.gz）
        mtx_files = list(dir_path.glob('matrix.mtx*'))
        if not mtx_files:
            raise FileNotFoundError(f"矩阵目录{dir_path}中未找到matrix.mtx文件")
        mtx_path = mtx_files[0]
        
        # 2. 读取稀疏矩阵（Matrix Market格式）
        try:
            mat = scipy.io.mmread(mtx_path).tocsc()  # 转为按列压缩（细胞为列）
        except Exception as e:
            raise RuntimeError(f"读取矩阵{mtx_path}失败: {str(e)}")
        matrices.append(mat)
        
        # 3. 查找基因列表文件（features.tsv或genes.tsv，支持.gz）
        feat_files = list(dir_path.glob('features.tsv*')) + list(dir_path.glob('genes.tsv*'))
        if not feat_files:
            raise FileNotFoundError(f"矩阵目录{dir_path}中未找到基因列表文件（features.tsv/genes.tsv）")
        feat_path = feat_files[0]
        
        # 4. 读取基因列表（第2列通常为基因名）
        compression = 'gzip' if str(feat_path).endswith('.gz') else None
        try:
            features = pd.read_csv(feat_path, sep='\t', header=None, compression=compression)
        except Exception as e:
            raise RuntimeError(f"读取基因列表{feat_path}失败: {str(e)}")
        current_genes = features[1].values  # 取第2列作为基因名
        
        # 5. 验证基因列表一致性（所有子库必须使用相同基因注释）
        if gene_list is None:
            gene_list = current_genes
        else:
            if not np.array_equal(gene_list, current_genes):
                raise ValueError(f"基因列表在{dir_path}与首个矩阵不一致，可能注释版本不同")
    
    # 合并所有矩阵（按列拼接，即合并所有细胞）
    combined_matrix = sparse.hstack(matrices)
    
    # 统计每个基因在所有细胞中的总表达量，计数表达量>0的基因数
    gene_total_expression = np.array(combined_matrix.sum(axis=1)).flatten()  # 按行求和（基因维度）
    total_genes = int((gene_total_expression > 0).sum())  # 表达量>0的基因数
    
    return total_genes


def main():
    # 解析命令行参数
    args = parse_args()
    
    # 1. 从splitcode summary文件提取总reads数（n_reads_max）
    print("步骤1/6：提取总测序reads数（n_reads_max）...")
    n_reads_max = extract_n_reads_max(args.inputfileSplitSummary)
    print(f"成功提取n_reads_max: {n_reads_max}")
    
    # 2. 读取Results_Summary.xls（子库基础指标）
    print("步骤2/6：读取子库统计结果（Results_Summary.xls）...")
    try:
        summary_df = pd.read_csv(args.inputfileSummary, sep='\t', index_col=0)
    except Exception as e:
        raise RuntimeError(f"读取{args.inputfileSummary}失败: {str(e)}")
    
    # 3. 计算库级汇总指标（基于Results_Summary和n_reads_max）
    print("步骤3/6：计算库级汇总指标...")
    library_metrics = summarize_library_metrics(summary_df, args.SplitSample, n_reads_max)
    
    # 4. 处理所有子库的CellReads.stats，生成细胞级汇总表
    print("步骤4/6：处理CellReads.stats，生成细胞级数据...")
    tmp_top_cell_paths = []  # 存储各子库处理后的临时文件路径
    for stat_path in args.inputfileCellReads:
        # 提取子库名称（从路径中解析，如star_outs_lib9 -> 9）
        lib_dir = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(stat_path))))
        lib_id = re.search(r'\d+$', lib_dir).group() if re.search(r'\d+$', lib_dir) else lib_dir
        
        # 获取该子库的估计细胞数
        if lib_id not in summary_df.columns:
            raise ValueError(f"子库{lib_id}在{args.inputfileSummary}中未找到对应列")
        est_cells = int(summary_df.loc['Estimated Number of Cells', lib_id])
        
        # 处理该子库的CellReads.stats
        print(f"  处理子库{lib_id}的CellReads.stats（估计细胞数：{est_cells}）...")
        tmp_path = process_cell_reads(stat_path, est_cells)
        tmp_top_cell_paths.append(tmp_path)
    
    # 合并所有子库的细胞级数据，输出至Cell_Summary.xls
    print("  合并所有细胞级数据...")
    cell_summary_dfs = [pd.read_csv(p, sep='\t', header=None) for p in tmp_top_cell_paths]
    combined_cell_df = pd.concat(cell_summary_dfs, ignore_index=True)
    combined_cell_df.to_csv(args.outputfile1, sep='\t', header=False, index=False)
    print(f"  细胞级汇总表已写入: {args.outputfile1}")
    
    # 清理细胞级处理的临时文件
    for tmp in tmp_top_cell_paths:
        if os.path.exists(tmp):
            os.remove(tmp)
    
    # 5. 计算细胞级指标（均值、中位数）和测序饱和度
    print("步骤5/6：计算细胞级指标和测序饱和度...")
    # 细胞级指标（第15列：GenicReads，第17列：UMI，第18列：Gene数）
    mean_genic_reads = int(round(combined_cell_df[14].mean()))  # 第15列（索引14）
    mean_umi = int(round(combined_cell_df[16].mean()))          # 第17列（索引16）
    mean_gene = int(round(combined_cell_df[17].mean()))         # 第18列（索引17）
    med_genic_reads = int(round(combined_cell_df[14].median()))
    med_umi = int(round(combined_cell_df[16].median()))
    med_gene = int(round(combined_cell_df[17].median()))
    
    # 测序饱和度 = (1 - 总唯一UMI数 / 总GenicReads数) * 100
    total_umi = combined_cell_df[16].sum()
    total_genic_reads = combined_cell_df[14].sum()
    sequencing_saturation = round((1 - total_umi / total_genic_reads) * 100, 2)
    
    # 6. 计算总基因数（Total Gene Detected）
    print("步骤6/6：合并10X矩阵，统计总基因数...")
    total_genes = compute_total_genes_covered(args.inputfileMatrixDir)
    print(f"  检测到的总基因数: {total_genes}")
    
    # 组装最终汇总指标（插入测序饱和度、细胞平均reads等）
    print("组装最终汇总指标...")
    metrics_list = list(library_metrics.items())  # 转换为列表以保持顺序
    
    # 插入测序饱和度（在有效Barcode指标后）
    valid_barcode_idx = next(i for i, (k, _) in enumerate(metrics_list) if k == 'Reads With Valid Barcodes')
    metrics_list.insert(valid_barcode_idx + 1, ('Sequencing Saturation', sequencing_saturation))
    
    # 插入细胞平均reads数（基于有效Barcode reads和总细胞数）
    mean_reads_valid = int(library_metrics['Reads With Valid Barcodes Num'] / library_metrics['Estimated Number of Cells'])
    frac_unique_idx = next(i for i, (k, _) in enumerate(metrics_list) if k == 'Fraction of Unique Reads in Cells')
    metrics_list.insert(frac_unique_idx + 1, ('Mean Reads per Cell (Valid Barcodes)', mean_reads_valid))
    
    # 插入细胞平均原始reads数（基于总reads数和总细胞数）
    mean_reads_raw = int(n_reads_max / library_metrics['Estimated Number of Cells'])
    metrics_list.insert(frac_unique_idx + 2, ('Mean Reads per Cell (Raw Reads)', mean_reads_raw))
    
    # 添加总基因数指标
    metrics_list.append(('Total Gene Detected', total_genes))
    
    # 生成最终汇总表并输出
    final_rows = [{'Sample': k, 'Merged': v} for k, v in metrics_list]
    final_rows.extend([
        {'Sample': 'Mean GenicReads per Cell', 'Merged': mean_genic_reads},
        {'Sample': 'Mean UMI per Cell', 'Merged': mean_umi},
        {'Sample': 'Mean Gene per Cell', 'Merged': mean_gene},
        {'Sample': 'Median GenicReads per Cell', 'Merged': med_genic_reads},
        {'Sample': 'Median UMI per Cell', 'Merged': med_umi},
        {'Sample': 'Median Gene per Cell', 'Merged': med_gene},
    ])
    final_df = pd.DataFrame(final_rows)
    final_df.to_csv(args.outputfile2, sep='\t', index=False)
    print(f"最终汇总指标已写入: {args.outputfile2}")
    
    print("所有分析完成！")


if __name__ == '__main__':
    main()