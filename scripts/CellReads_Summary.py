# -*- coding: utf-8 -*-

import csv
import os
import sys

# 检查命令行参数
if len(sys.argv) != 4:
    print("Usage: python CellReads_Summary.py <input_file> <output_file> <header_value>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]
header_value = sys.argv[3]

# 第一步：去除第二行，生成临时文件
with open(input_file, 'r') as f:
    lines = f.readlines()

# 保留标题行和第三行及之后的内容
filtered_lines = [lines[0]] + lines[2:]

# 写入临时文件
temp_file = "CellReads.stats.tmp"
with open(temp_file, 'w') as f:
    f.writelines(filtered_lines)

# 第二步：计算各列求和并生成结果表
# 初始化各列求和变量
genomeU_sum = 0
genomeM_sum = 0
exonic_sum = 0
intronic_sum = 0
exonicAS_sum = 0
intronicAS_sum = 0

# 读取临时文件并计算列总和
with open(temp_file, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        genomeU_sum += int(row['genomeU'])
        genomeM_sum += int(row['genomeM'])
        exonic_sum += int(row['exonic'])
        intronic_sum += int(row['intronic'])
        exonicAS_sum += int(row['exonicAS'])
        intronicAS_sum += int(row['intronicAS'])

# 计算衍生指标
antisense = exonicAS_sum + intronicAS_sum
intergenic = genomeU_sum + genomeM_sum - exonic_sum - intronic_sum - antisense

# 计算分母总和
total = intergenic + antisense + exonic_sum + intronic_sum

# 计算各比例（处理除零错误）
if total != 0:
    exonic_frac = round(exonic_sum / total * 100, 2)
    intronic_frac = round(intronic_sum / total * 100, 2)
    intergenic_frac = round(intergenic / total * 100, 2)
    antisense_frac = round(antisense / total * 100, 2)
else:
    exonic_frac = intronic_frac = intergenic_frac = antisense_frac = 0.00

# 构建结果数据
results = [
    ('Uniquely Mapped Reads', genomeU_sum),
    ('Multi-Mapped Reads', genomeM_sum),
    ('Mapped Reads Assigned To Exonic Regions', exonic_sum),
    ('Mapped Reads Assigned To Intronic Regions', intronic_sum),
    ('Mapped Reads Assigned To Intergenic Regions', intergenic),
    ('Mapped Reads Assigned Antisense To Gene', antisense),
    ('Mapped Reads Assigned To Exonic Regions Fraction', exonic_frac),
    ('Mapped Reads Assigned To Intronic Regions Fraction', intronic_frac),
    ('Mapped Reads Assigned To Intergenic Regions Fraction', intergenic_frac),
    ('Mapped Reads Assigned Antisense To Gene Fraction', antisense_frac),
]

# 写入结果文件，header_value 由命令行参数提供
with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['Statistic', header_value])
    for row in results:
        writer.writerow(row)

# 删除临时文件
os.remove(temp_file)

print(f"处理完成，结果已保存至 {output_file}")
