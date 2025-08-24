import argparse
import pandas as pd
import numpy as np
import scanpy as sc
import re
from pathlib import Path
import os

# ------------------ 函数定义 ------------------
def parse_index_spec(spec, total_possible_indices=None):
    """
    解析索引规范（兼容三种格式，0-based索引）：
    1. 等差数列：[start,end,step] 如[0,88,8]（0-based，包含起止，步长8）
    2. 范围格式：0-4（0-based，包含起止）
    3. 单个索引：2（0-based）
    
    返回：0-based索引列表（内部直接使用）
    """
    indices = []
    
    # 处理等差数列格式 [start,end,step]
    if spec.startswith('[') and spec.endswith(']'):
        try:
            content = spec[1:-1].strip()
            start, end, step = map(int, content.split(','))
            if step <= 0:
                raise ValueError(f"步长必须>0，当前值：{step}")
            if start > end:
                start, end = end, start
            # 生成0-based索引
            current = start
            while current <= end:
                indices.append(current)
                current += step
        except Exception as e:
            raise ValueError(f"无效等差数列格式 '{spec}': {e}")
        return indices
    
    # 处理范围或单个索引格式
    parts = spec.split(',')
    for part in parts:
        part = part.strip()
        if '-' in part:
            # 范围格式：0-4 → 0-based（包含起止）
            try:
                a, b = map(int, part.split('-'))
                if a < 0 or b < a:
                    raise ValueError(f"范围 [{a}-{b}] 无效（0-based索引需满足a≥0且b≥a）")
                indices.extend(range(a, b + 1))  # 包含结束位置
            except:
                raise ValueError(f"无效范围格式 '{part}'，应为0-based的start-end（如0-4）")
        else:
            # 单个索引：2 → 0-based
            try:
                idx = int(part)
                if idx < 0:
                    raise ValueError(f"索引 {idx} 无效（0-based索引需≥0）")
                indices.append(idx)
            except:
                raise ValueError(f"无效单个索引 '{part}'，应为0-based整数（如2）")
    
    # 去重并保持顺序
    seen = set()
    return [x for x in indices if x not in seen and not seen.add(x)]

# ------------------ 主流程 ------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='SampleID split and summary（基于mapping.txt的bc2-xxx索引分配样本）',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('--inputfileCellReads', nargs='+', help='CellReads.stats 文件列表', required=True)
    parser.add_argument('--inputfileBarcodelist', help='Barcode映射文件（mapping.txt，tab分隔）', required=True)
    parser.add_argument('--inputfileCellSummaryFilter', help='过滤后带SampleID的CellSummary文件', required=True)
    parser.add_argument('--splitCB', nargs='+', required=True,
                        metavar='SPEC',
                        help='0-based索引规范列表，支持：\n'
                             '1. 等差数列：[起始,结束,步长] 如[0,88,8]\n'
                             '2. 范围格式：起始-结束 如0-4\n'
                             '3. 单个索引：2\n'
                             '参数格式：双引号包裹，空格分隔多个规范，如 "[0,88,8] [1,89,8]"')
    parser.add_argument('--splitSample', nargs='+', required=True,
                        metavar='SAMPLE',
                        help='与splitCB数量对应的样本名称列表（顺序严格匹配）\n'
                             '参数格式：双引号包裹，空格分隔多个样本，如 "A A1 B"')
    parser.add_argument('--inputfileCellMatrix', nargs='+', help='filtered_feature_bc_matrix目录列表', required=True)
    parser.add_argument('--outputfileCellSummaryRaw', help='输出带SampleID的原始表', required=True)
    parser.add_argument('--outputfileSampleSummary', help='输出样本汇总表', required=True)
    args = parser.parse_args()

    # ====================== 参数校验 ======================
    if len(args.splitCB) != len(args.splitSample):
        raise ValueError(f"splitCB({len(args.splitCB)})与splitSample({len(args.splitSample)})数量不匹配")

    # ====================== 路径处理 ======================
    MATRIX_BASE_DIR = Path("07-CB2_Matrix")
    MATRIX_BASE_DIR.mkdir(parents=True, exist_ok=True)
    
    processed_matrix_paths = []
    for matrix_name in args.inputfileCellMatrix:
        path_obj = Path(matrix_name)
        if path_obj.is_absolute():
            processed_path = path_obj
        else:
            processed_path = MATRIX_BASE_DIR / matrix_name
        if not processed_path.exists():
            raise FileNotFoundError(f"矩阵目录不存在: {processed_path}")
        processed_matrix_paths.append(str(processed_path))
    args.inputfileCellMatrix = processed_matrix_paths

    # ====================== 读取并解析mapping.txt ======================
    # 从mapping.txt提取Barcode及其对应的索引（从第二列的bc2-xxx解析）
    barcode_to_index = {}  # key: 16bp Barcode (第一列), value: 0-based索引 (从bc2-xxx提取)
    index_pattern = re.compile(r'bc2-(\d+)')  # 匹配"bc2-数字"的正则表达式

    with open(args.inputfileBarcodelist) as f:
        line_num = 0
        for line in f:
            line_num += 1
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')  # tab分隔
            if len(parts) < 2:
                raise ValueError(f"mapping.txt第{line_num}行格式错误，至少需要两列（tab分隔）")
            
            barcode = parts[0].strip()  # 第一列：16bp Barcode（完整使用，不截取）
            second_col = parts[1].strip()  # 第二列：包含bc2-xxx的信息
            
            # 从第二列提取bc2-后的数字作为索引
            match = index_pattern.search(second_col)
            if not match:
                raise ValueError(f"mapping.txt第{line_num}行第二列未找到bc2-xxx格式（如bc2-0），内容：{second_col}")
            
            try:
                idx = int(match.group(1))  # 提取数字部分并转换为整数
                if idx < 0:
                    raise ValueError(f"mapping.txt第{line_num}行索引为负数：{idx}")
            except ValueError:
                raise ValueError(f"mapping.txt第{line_num}行bc2-后不是有效整数：{match.group(1)}")
            
            if barcode in barcode_to_index:
                raise ValueError(f"mapping.txt中Barcode '{barcode}' 重复出现（第{line_num}行）")
            barcode_to_index[barcode] = idx
            # 调试信息：打印前5条解析结果
            if line_num <= 5:
                print(f"解析mapping.txt第{line_num}行：Barcode={barcode} → 索引={idx}（来自{second_col}）")

    if not barcode_to_index:
        raise ValueError(f"mapping.txt中未解析到有效Barcode数据：{args.inputfileBarcodelist}")
    print(f"成功解析mapping.txt，共获取 {len(barcode_to_index)} 个Barcode及其索引")

    # ====================== 解析splitCB规范为样本-索引范围映射 ======================
    # 建立索引到样本的映射：index_to_sample[idx] = 样本名
    index_to_sample = {}
    for spec_idx, (spec, sample) in enumerate(zip(args.splitCB, args.splitSample)):
        # 解析当前spec对应的索引列表（0-based）
        indices = parse_index_spec(spec)
        # 检查这些索引是否已被分配给其他样本
        for idx in indices:
            if idx in index_to_sample:
                raise ValueError(f"索引 {idx} 被分配到多个样本：{index_to_sample[idx]} 和 {sample}（来自规范 {spec}）")
            index_to_sample[idx] = sample
        print(f"样本 {sample} 对应的索引范围（来自规范 {spec}）：{indices[:5]}...（共{len(indices)}个）")

    # ====================== 建立Barcode到样本的映射 ======================
    barcode_to_sample = {}
    for barcode, idx in barcode_to_index.items():
        # 根据Barcode的索引查找对应的样本
        sample = index_to_sample.get(idx, 'Unknown')
        barcode_to_sample[barcode] = sample
        # 调试信息：打印前5条映射结果
        if len(barcode_to_sample) <= 5:
            print(f"Barcode={barcode}（索引={idx}）→ 样本={sample}")

    # ====================== 处理CellReads.stats文件 ======================
    # 合并所有CellReads.stats文件
    df_list = []
    for fn in args.inputfileCellReads:
        df = pd.read_csv(fn, sep='\t', header=0, skiprows=[1], dtype={'CB': str})
        df_list.append(df)
    df_all = pd.concat(df_list, axis=0, ignore_index=True)
    df_all.sort_values(by=df_all.columns[16], ascending=False, inplace=True)

    # 直接使用CB列的完整Barcode分配SampleID（无需截取）
    def assign_sample(cb):
        """根据完整Barcode查找对应的样本"""
        return barcode_to_sample.get(cb, 'Unknown')  # 不截取，直接用完整CB匹配

    df_all['SampleID'] = df_all['CB'].apply(assign_sample)
    df_all.to_csv(args.outputfileCellSummaryRaw, sep='\t', index=False)
    print(f"已生成带SampleID的原始表: {args.outputfileCellSummaryRaw}")

    # ====================== 计算各项指标 ======================
    # 基于CellReads原始表计算整体mapping指标
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

    # 处理过滤后CellSummary表
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
    
    # 计算测序饱和度（除零保护）
    seq_sat = pd.Series(dtype=str)
    for sample in args.splitSample:
        valid_count = countedU_sum_filt.get(sample, 0)
        seq_sat[sample] = f"{(1 - nUMIunique_sum_filt.get(sample, 0)/valid_count)*100:.2f}" if valid_count > 0 else "0.00"

    # 计算基因覆盖数
    genes_covered = {}
    for matrix_dir, sample in zip(args.inputfileCellMatrix, args.splitSample):
        try:
            print(f"处理样本 {sample} 的矩阵: {matrix_dir}")
            adata = sc.read_10x_mtx(matrix_dir, cache=False)
            gene_counts = np.array(adata.X.sum(axis=0)).flatten()
            genes_covered[sample] = int((gene_counts > 0).sum())
            print(f"样本 {sample}: 检测到 {genes_covered[sample]} 个基因")
        except Exception as e:
            raise ValueError(f"读取矩阵文件 {matrix_dir} 失败: {e}")

    # 计算平均reads数（除零保护）
    lib_mean = pd.Series(dtype=int)
    for sample in args.splitSample:
        rv = reads_valid.get(sample, 0)
        ec = est_cells.get(sample, 0)
        lib_mean[sample] = round(rv / ec) if ec > 0 else 0

    # ====================== 整合并输出结果 ======================
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

    summary_data = {sample: {name: vals.get(sample, np.nan) if isinstance(vals, pd.Series) else vals for name, vals in metrics} for sample in args.splitSample}
    summary_df = pd.DataFrame(summary_data).T
    summary_df.to_csv(args.outputfileSampleSummary, sep='\t', index=True, header=True)
    print(f"已生成样本汇总表: {args.outputfileSampleSummary}")