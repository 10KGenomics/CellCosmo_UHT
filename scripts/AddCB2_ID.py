# -*- coding: utf-8 -*-

"""
AddCB2_ID.py

运行说明：
1. 读取多个 CellReads.stats 文件；
2. 去除第二行 (CBnotInPasslist)；
3. 按 nUMIunique 降序排序，提取前 4500 个 cell（加表头共4501行）；
4. 从 CB 列中提取第二段 barcode，与大标签拼接生成 CB2_ID；
5. 合并所有文件为一个总表 Cell_AddCB2ID_XXX_Summary.xls；
6. 按 CB2_ID 分组，统计若干列的求和、中位数、平均数、百分比以及 cell 数；
7. 保存统计结果为 CB2ID_XXX_Summary.xls；
"""

import pandas as pd
import os
import argparse
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description='Add CB2_ID column and summarize statistics by sample.')
    parser.add_argument('--inputfile', nargs='+', required=True, help='List of CellReads.stats files')
    parser.add_argument('--outputfile1', required=True, help='Output file with CB2_ID added and merged')
    parser.add_argument('--group_col', default='CB2_ID', help='Grouping column, default is CB2_ID')
    parser.add_argument('--stat_cols', nargs='+', required=True, help='Columns to summarize')
    parser.add_argument('--outputfile2', required=True, help='Summary statistics file by CB2_ID group')
    return parser.parse_args()

def extract_tag_from_path(path):
    # 提取大标签，如从路径 02-STARsolo/2450_1_XXX-Solo.out/... 中提取 2450_1
    basename = os.path.basename(os.path.dirname(os.path.dirname(path)))
    return '_'.join(basename.split('_')[:1])

def add_CB2_ID_column(df, tag):
    df = df.copy()
    df['CB2_ID'] = df['CB'].apply(lambda x: f"{tag}_{x.split('_')[1]}" if '_' in x else 'NA')
    return df

def read_and_process_file(filepath):
    tag = extract_tag_from_path(filepath)
    with open(filepath) as f:
        lines = f.readlines()
    # 去掉第二行（CBnotInPasslist）
    header = lines[0].strip().split('\t')
    lines_filtered = [lines[0]] + lines[2:]

    # 读入 DataFrame
    df = pd.DataFrame([x.strip().split('\t') for x in lines_filtered[1:]], columns=header)
    
    # 转换数值列为 int
    for col in df.columns[1:]:
        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)

    # 按 nUMIunique 排序，取前 4500 个 cell（加 header 共 1001 行）
    df = df.sort_values(by='nUMIunique', ascending=False).head(4500)

    # 加 CB2_ID 列
    df = add_CB2_ID_column(df, tag)
    return df

def main():
    args = parse_args()
    all_dfs = []

    for file in args.inputfile:
        df = read_and_process_file(file)
        all_dfs.append(df)

    combined_df = pd.concat(all_dfs, ignore_index=True)
    combined_df.to_csv(args.outputfile1, sep='\t', index=False)

    # 按 CB2_ID 分组统计
    group_stats = []
    group_col = args.group_col
    stat_cols = args.stat_cols

    grouped = combined_df.groupby(group_col)
    for group, gdf in grouped:
        group_result = {
            group_col: group,
            'CellNum': len(gdf)
        }
        total_sum = gdf[stat_cols].sum()
        total_median = gdf[stat_cols].median()
        total_mean = gdf[stat_cols].mean()

        for col in stat_cols:
            group_result[f"{col}_sum"] = int(total_sum[col])
            group_result[f"{col}_median"] = int(total_median[col])
            group_result[f"{col}_mean"] = int(total_mean[col])

        group_stats.append(group_result)

    result_df = pd.DataFrame(group_stats)

    # 计算百分比列
    for col in stat_cols:
        total = result_df[f"{col}_sum"].sum()
        result_df[f"{col}_percent"] = (result_df[f"{col}_sum"] / total * 100).round(2)

    result_df.to_csv(args.outputfile2, sep='\t', index=False)

if __name__ == '__main__':
    main()
