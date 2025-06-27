# -*- coding: utf-8 -*-
"""
Combined_Sample.py

将多个样本 Summary 文件合并为一个表格，输出文件格式为制表符分隔的文本文件。

用法:
    python Combined_Sample.py <file1> <file2> ... <Output_File>

示例:
    python Combined_Sample.py \
      02-STARsolo/16to1_GATCGGTA-Solo.out/GeneFull_Ex50pAS/16to1_GATCGGTA-Summary.xls \
      02-STARsolo/16to1_GGACGTAT-Solo.out/GeneFull_Ex50pAS/16to1_GGACGTAT-Summary.xls \
      Results_Summary.xls

功能说明：
    1. 原脚本中直接读取 Summary 文件中已有的 "Mean Reads per Cell" 和 "Median Reads per Cell"，
       实际对应的是基因组映射后的基因 Reads（Genic Reads）统计。为避免误解，
       将其重命名为：
           "Mean GenicReads per Cell" 和 "Median GenicReads per Cell"。
    2. 新增对每个样本的整体 Mean Reads per Cell 计算：
       使用 Reads With Valid Barcodes Num 除以 Estimated Number of Cells，并取整，
       命名为 "Mean Reads per Cell"。不计算中位数。
    3. 输出表中，将新计算的 Mean Reads per Cell 插入在
       "Fraction of Unique Reads in Cells" 之后；
       接着放置重命名的 Mean/Median GenicReads per Cell。
    4. 详细注释已在代码对应部分标注。
"""

import csv
import sys

if len(sys.argv) < 3:
    print("Usage: python Combined_Sample.py <file1> <file2> ... <Output_File>")
    sys.exit(1)

# 输入文件列表，最后一项为输出文件名
input_files = sys.argv[1:-1]
output_file = sys.argv[-1]

# 样本名称列表（按输入顺序）和样本数据字典
sample_names = []      # 存放样本名
sample_data = {}       # {sample_name: {statistic: value, ...}}

# 基准统计项顺序，来自第一个文件
order_stats = []

# 读取每个 Summary 文件
for file in input_files:
    with open(file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        if len(header) < 2:
            print(f"文件 {file} 的格式不正确，缺少样本名")
            sys.exit(1)
        sample_name = header[1].strip()
        sample_names.append(sample_name)
        current_dict = {}
        for row in reader:
            if len(row) < 2:
                continue
            stat = row[0].strip()
            value = row[1].strip()
            current_dict[stat] = value
        sample_data[sample_name] = current_dict
        if not order_stats:
            order_stats = list(current_dict.keys())

# 收集所有统计项（包括额外项）
all_stats = set(order_stats)
for s in sample_names:
    all_stats.update(sample_data[s].keys())
# 保留基准顺序，其它按字母排序追加
extra_stats = sorted(list(all_stats - set(order_stats)))
combined_stats = order_stats + extra_stats

# 重命名和插入说明：
orig_mean = "Mean GenicReads per Cell"               # 原文件中的 Mean Reads per Cell (实际上为 GenicReads)   #二次优化，上一个脚本已更改
orig_median = "Median GenicReads per Cell"           # 原文件中的 Median Reads per Cell (实际上为 GenicReads) #二次优化，上一个脚本已更改
new_mean = "Mean Reads per Cell(based on Reads With Valid Barcodes Num)"                     # 新增：整体 Mean Reads per Cell
genic_mean = "Mean GenicReads per Cell"              # 重命名后：Genic Reads    均值
genic_median = "Median GenicReads per Cell"          # 重命名后：Genic Reads    中位数

# 构造新的统计项顺序：跳过原始 Mean/Median，
# 在 "Fraction of Unique Reads in Cells" 后插入 new_mean、genic_mean、genic_median
new_combined_stats = []
for stat in combined_stats:
    # 跳过原始名
    if stat in (orig_mean, orig_median):
        continue
    new_combined_stats.append(stat)
    if stat == "Fraction of Unique Reads in Cells":
        # 插入新计算的整体 Mean Reads per Cell
        new_combined_stats.append(new_mean)
        # 插入重命名的基因 Reads 统计
        new_combined_stats.append(genic_mean)
        new_combined_stats.append(genic_median)

# 构造输出表格
output_table = []
header_row = ["Statistic"] + sample_names
output_table.append(header_row)

for stat in new_combined_stats:
    row = [stat]
    for sname in sample_names:
        if stat == new_mean:
            # 计算：Reads With Valid Barcodes Num / Estimated Number of Cells，保留整数
            num = sample_data[sname].get("Reads With Valid Barcodes Num", "0")
            cells = sample_data[sname].get("Estimated Number of Cells", "1")
            try:
                value = str(int(int(num) / int(cells)))
            except Exception:
                value = ""
        elif stat == genic_mean:
            # 使用原 Mean Reads per Cell 值
            value = sample_data[sname].get(orig_mean, "")
        elif stat == genic_median:
            # 使用原 Median Reads per Cell 值
            value = sample_data[sname].get(orig_median, "")
        else:
            # 其他统计项按原始读取
            value = sample_data[sname].get(stat, "")
        row.append(value)
    output_table.append(row)

# 写出结果
with open(output_file, 'w', newline='') as fout:
    writer = csv.writer(fout, delimiter='\t')
    writer.writerows(output_table)

print(f"合并完成，结果已保存至 {output_file}")
