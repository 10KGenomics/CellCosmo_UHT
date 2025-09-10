# -*- coding: utf-8 -*-
import csv
import sys

# 检查命令行参数
if len(sys.argv) != 5:
    print("Usage: python Summary_merge.py <CellReads_Summary.xls> <Summary.csv> <Output_file> <Sample_name>")
    sys.exit(1)

cellreads_file = sys.argv[1]   # e.g., "16to1_GATCGGTA-CellReads_Summary.xls"
csv_file       = sys.argv[2]   # e.g., "Summary.csv"
output_file    = sys.argv[3]   # e.g., "16to1_GATCGGTA-Summary.xls"
sample_name    = sys.argv[4]   # e.g., "16to1_GATCGGTA"

# ----------------------------
# 读取 CSV 文件（Summary.csv）
# CSV 文件格式：每行格式为 Key,Value
csv_data = {}
with open(csv_file, 'r') as f:
    reader = csv.reader(f, delimiter=',')
    for row in reader:
        if len(row) >= 2:
            key = row[0].strip()
            value = row[1].strip()
            csv_data[key] = value

# 转换并提取CSV中所需数据（注意：部分值为小数，需要转为float）
num_reads = int(csv_data["Number of Reads"])
rwb_fraction = float(csv_data["Reads With Valid Barcodes"])  
rwb_pct = round(rwb_fraction * 100, 2)
rwb_num = int(round(num_reads * rwb_fraction))
sequencing_saturation_pct = round(float(csv_data["Sequencing Saturation"]) * 100, 2)
# q30_cb_umi_pct = round(float(csv_data["Q30 Bases in CB+UMI"]) * 100, 2)
q30_rna_pct = round(float(csv_data["Q30 Bases in RNA read"]) * 100, 2)
frac_unique_reads_in_cells_pct = round(float(csv_data["Fraction of Unique Reads in Cells"]) * 100, 2)
estimated_cells = int(csv_data["Estimated Number of Cells"])
unique_reads_in_cells = int(csv_data["Unique Reads in Cells Mapped to GeneFull_Ex50pAS"])
mean_reads_per_cell = csv_data["Mean Reads per Cell"]
median_reads_per_cell = csv_data["Median Reads per Cell"]
umis_in_cells = csv_data["UMIs in Cells"]
mean_umi_per_cell = csv_data["Mean UMI per Cell"]
median_umi_per_cell = csv_data["Median UMI per Cell"]
mean_genes_per_cell = csv_data["Mean GeneFull_Ex50pAS per Cell"]
median_genes_per_cell = csv_data["Median GeneFull_Ex50pAS per Cell"]
total_genes_detected = csv_data["Total GeneFull_Ex50pAS Detected"]

# ----------------------------
# 读取 CellReads_Summary 文件
cellreads_data = {}
with open(cellreads_file, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        key = row["Statistic"].strip()
        cellreads_data[key] = row[next(iter(row.keys() - {"Statistic"}))].strip()

uniquely_mapped = int(cellreads_data["Uniquely Mapped Reads"])
multi_mapped = int(cellreads_data["Multi-Mapped Reads"])
mapped_exonic = int(cellreads_data["Mapped Reads Assigned To Exonic Regions"])
mapped_intronic = int(cellreads_data["Mapped Reads Assigned To Intronic Regions"])
mapped_intergenic = int(cellreads_data["Mapped Reads Assigned To Intergenic Regions"])
mapped_antisense = int(cellreads_data["Mapped Reads Assigned Antisense To Gene"])
exonic_frac = float(cellreads_data["Mapped Reads Assigned To Exonic Regions Fraction"])
intronic_frac = float(cellreads_data["Mapped Reads Assigned To Intronic Regions Fraction"])
intergenic_frac = float(cellreads_data["Mapped Reads Assigned To Intergenic Regions Fraction"])
antisense_frac = float(cellreads_data["Mapped Reads Assigned Antisense To Gene Fraction"])

# ----------------------------
uniquely_mapped_fraction = round((uniquely_mapped / rwb_num) * 100, 2)
multi_mapped_fraction = round((multi_mapped / rwb_num) * 100, 2)

# ----------------------------
results = [
    ("Number of Reads", num_reads),
    # ("Reads With Valid Barcodes", f"{rwb_pct}"),
    ("Reads With Valid Barcodes Num", rwb_num),
    ("Sequencing Saturation", f"{sequencing_saturation_pct}"),
    # ("Q30 Bases in CB+UMI", f"{q30_cb_umi_pct}"),
    ("Q30 Bases in RNA read", f"{q30_rna_pct}"),
    ("Uniquely Mapped Reads", uniquely_mapped),
    ("Uniquely Mapped Reads Fraction", f"{uniquely_mapped_fraction}"),
    ("Multi-Mapped Reads", multi_mapped),
    ("Multi-Mapped Reads Fraction", f"{multi_mapped_fraction}"),
    ("Mapped Reads Assigned To Exonic Regions", mapped_exonic),
    ("Mapped Reads Assigned To Intronic Regions", mapped_intronic),
    ("Mapped Reads Assigned To Intergenic Regions", mapped_intergenic),
    ("Mapped Reads Assigned Antisense To Gene", mapped_antisense),
    ("Mapped Reads Assigned To Exonic Regions Fraction", exonic_frac),
    ("Mapped Reads Assigned To Intronic Regions Fraction", intronic_frac),
    ("Mapped Reads Assigned To Intergenic Regions Fraction", intergenic_frac),
    ("Mapped Reads Assigned Antisense To Gene Fraction", antisense_frac),
    ("Estimated Number of Cells", estimated_cells),
    ("Unique Reads in Cells Mapped to Gene", unique_reads_in_cells),
    ("Fraction of Unique Reads in Cells", f"{frac_unique_reads_in_cells_pct}"),
    ("Mean GenicReads per Cell", mean_reads_per_cell),
    ("Median GenicReads per Cell", median_reads_per_cell),
    ("UMIs in Cells", umis_in_cells),
    ("Mean UMI per Cell", mean_umi_per_cell),
    ("Median UMI per Cell", median_umi_per_cell),
    ("Mean Gene per Cell", mean_genes_per_cell),
    ("Median Gene per Cell", median_genes_per_cell),
    ("Total Gene Detected", total_genes_detected)
]

# ----------------------------
with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(["Statistic", sample_name])
    for item in results:
        writer.writerow(item)

print(f"合并完成，结果已保存至 {output_file}")
