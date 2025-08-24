import glob
import os
import base64
from jinja2 import Template
import pandas as pd
import re
from natsort import natsorted  # 导入自然排序库

# ---------------------------
# 1. 统计项 & 图表指标定义（完整保留）
# ---------------------------
stat_descriptions = {
    "Number of Reads": "原始测序下机reads数",
    "Reads With Valid Barcodes": "有效Barcode的reads占比",
    "Reads With Valid Barcodes Num": "有效Barcode的reads数（有效文库reads数）",
    "Sequencing Saturation": "测序饱和度",
    "Q30 Bases in CB+UMI": "CB+UMI序列区域的Q30碱基占比",
    "Q30 Bases in RNA read": "RNA read序列区域的Q30碱基占比",
    "Uniquely Mapped Reads": "唯一比对的reads数",
    "Uniquely Mapped Reads Fraction": "唯一比对的reads占比",
    "Multi-Mapped Reads": "多重比对的reads数",
    "Multi-Mapped Reads Fraction": "多重比对的reads占比",
    "Mapped Reads Assigned To Exonic Regions": "外显子区域映射的reads数",
    "Mapped Reads Assigned To Intronic Regions": "内含子区域映射的reads数",
    "Mapped Reads Assigned To Intergenic Regions": "基因间区映射的reads数",
    "Mapped Reads Assigned Antisense To Gene": "基因反义链区域映射的reads数",
    "Mapped Reads Assigned To Exonic Regions Fraction": "外显子区域映射reads占比",
    "Mapped Reads Assigned To Intronic Regions Fraction": "内含子区域映射reads占比",
    "Mapped Reads Assigned To Intergenic Regions Fraction": "基因间区映射reads占比",
    "Mapped Reads Assigned Antisense To Gene Fraction": "基因反义链区域映射reads占比",
    "Estimated Number of Cells": "预估细胞数",
    "Unique Reads in Cells Mapped to Gene": "细胞中映射到基因的唯一reads数",
    "Fraction of Unique Reads in Cells": "细胞中唯一比对reads的占比",
    "Mean Reads per Cell(based on Number of RawReads)": "平均每个细胞的原始测序reads数分布值",
    "Mean Reads per Cell(based on Reads With Valid Barcodes Num)": "平均每个细胞的有效文库reads数分布值",
    "Mean GenicReads per Cell": "每个细胞的平均基因reads数",
    "Median GenicReads per Cell": "每个细胞的中位基因reads数",
    "UMIs in Cells": "细胞中的UMI总数",
    "Mean UMI per Cell": "每个细胞的平均UMI数",
    "Median UMI per Cell": "每个细胞的中位UMI数",
    "Mean Gene per Cell": "每个细胞的平均基因数",
    "Median Gene per Cell": "每个细胞的中位基因数",
    "Total Gene Detected": "检测到的总基因数"
}

CARD_METRICS = [
    "Estimated Number of Cells",
    "Median UMI per Cell",
    "Median Gene per Cell",
    "Total Gene Detected"
]

PIE_METRICS = [
    "Mapped Reads Assigned To Exonic Regions Fraction",
    "Mapped Reads Assigned To Intronic Regions Fraction",
    "Mapped Reads Assigned To Intergenic Regions Fraction",
    "Mapped Reads Assigned Antisense To Gene Fraction"
]

# ---------------------------
# 2. Logo读取（完整保留）
# ---------------------------
def get_logo_base64():
    try:
        with open("/home/rs1/3-script/10k_scRNA_UHT/10K-logo.png", "rb") as f:
            return base64.b64encode(f.read()).decode("utf-8")
    except FileNotFoundError:
        print("错误：请将 10K-logo.png 放在脚本同一目录！")
        exit(1)

logo_base64 = get_logo_base64()

# ---------------------------
# 3. 数据解析（核心修正区，此处保留原有逻辑，仅优化前端展示）
# ---------------------------
data_dict = {}                # Library单样本统计: {sample: [(stat, val), ...]}
barcode_rank_data = {}        # Library单样本Barcode: {sample: [{"x":1, "y":1000}, ...]}
sample_summary_data = {}      # Sample统计: {sample: [(stat, val), ...]}
sample_barcode_rank_data = {} # Sample Barcode: {sample: [{"x":1, "y":1000}, ...]}
merge_stats = []              # Merge统计: [(stat, val), ...]
merge_card_data = {}          # Merge卡片数据
merge_pie_data = {}           # Merge饼图数据
merge_barcode_rank_data = []  # Merge Barcode数据

# ---------------------------
# 3.1 解析Library单样本统计（03-CellReadsSummary/*-Summary.xls）
# 修正：统一按Tab分割，严格处理有效行（保留原有逻辑）
# ---------------------------
# 获取所有统计文件
stats_files = glob.glob("03-CellReadsSummary/*-Summary.xls")
print(f"发现 {len(stats_files)} 个Library单样本统计文件")

# 使用自然排序对文件列表进行排序
stats_files = natsorted(stats_files)

# 存储所有样本数据
data_dict = {}

# 处理每个文件
for file_path in stats_files:
    try:
        # 提取样本名称（在排序后提取）
        file_name = os.path.basename(file_path)
        sample_name = file_name.replace("-Summary.xls", "")
        
        # 读取文件内容
        with open(file_path, 'r') as f:
            lines = [line.rstrip('\n') for line in f.readlines()]
        
        # 筛选有效行
        valid_lines = []
        for line in lines:
            stripped = line.strip()
            # 保留标题行和有效数据行（包含制表符分隔的两列）
            if stripped.startswith("Statistic") or (len(stripped.split('\t')) == 2):
                valid_lines.append(line)
        
        if not valid_lines:
            print(f"警告：{file_path} 无有效数据，跳过")
            continue
        
        # 提取统计指标
        stats = []
        for line in valid_lines[1:]:  # 跳过标题行
            stripped = line.strip()
            if not stripped:
                continue
            stat_parts = stripped.split('\t', 1)
            if len(stat_parts) != 2:
                continue
            stat_name = stat_parts[0].strip()
            stat_value = stat_parts[1].strip()
            stats.append((stat_name, stat_value))
        
        # 保存到字典
        data_dict[sample_name] = stats
        print(f"✓ 解析Library样本 {sample_name}（{len(stats)} 项）")
    
    except Exception as e:
        print(f"错误：处理 {file_path} 失败：{str(e)}")

# 现在data_dict中的样本已经按名称自然排序

# ---------------------------
# 3.2 解析Library单样本Barcode（CellReads.stats.tmp2）
# 修正：添加数据抽样逻辑（保留原有逻辑，仅优化前端数据量）
# ---------------------------
SAMPLING_INTERVAL = 50  # 抽样间隔，每5个点取一个

for sample in data_dict.keys():
    # 使用原始的CellReads.stats文件路径
    barcode_file = f"02-STARsolo/star_outs_lib{sample}/Solo.out/GeneFull_Ex50pAS/CellReads.stats"
    try:
        with open(barcode_file, 'r') as f:
            # 跳过前两行：第一行是标题，第二行是"CBnotInPasslist"
            next(f)  # 跳过第一行标题
            next(f)  # 跳过第二行CBnotInPasslist
            lines = f.readlines()
        
        umi_counts = []
        for line in lines:
            parts = line.strip().split()
            # 确保有足够的列（至少17列）
            if len(parts) >= 17:
                try:
                    # 第17列（索引16）是countedU列
                    umi = int(parts[16])
                    umi_counts.append(umi)
                except ValueError:
                    # 处理可能的非整数值
                    print(f"警告：{sample} Barcode文件某行第17列非整数，跳过")
                    continue
        
        if not umi_counts:
            print(f"警告：{sample} 无有效Barcode数据，图表为空")
            barcode_rank_data[sample] = []
            continue
            
        # 按UMI计数降序排序
        umi_counts.sort(reverse=True)
        
        # 采样数据点（减少绘图数据量）
        sampled_umi_counts = umi_counts[::SAMPLING_INTERVAL]
        
        # 生成排名数据（x=排名，y=UMI计数）
        rank_data = [{"x": i*SAMPLING_INTERVAL+1, "y": val} for i, val in enumerate(sampled_umi_counts)]
        barcode_rank_data[sample] = rank_data
        
        print(f"✓ 解析 {sample} Barcode数据（原始{len(umi_counts)}点，抽样后{len(rank_data)}点）")
    
    except FileNotFoundError:
        print(f"错误：{barcode_file} 不存在，无法绘制Barcode图")
        barcode_rank_data[sample] = []
    except Exception as e:
        print(f"错误：处理 {sample} Barcode数据失败：{str(e)}")
        barcode_rank_data[sample] = []

# ---------------------------
# 3.3 解析Merge统计数据（05-Combined/2-Cell3_Merge_Summary.xls，Tab分隔）
# 修正：纯文本解析（保留原有逻辑）
# ---------------------------
merge_stats_file = "05-Combined/2-Cell3_Merge_Summary.xls"

def parse_merge_summary(file_path):
    if not os.path.exists(file_path):
        print(f"错误：Merge数据文件 {file_path} 不存在！")
        return []
    
    try:
        with open(file_path, 'r') as f:
            lines = [line.rstrip('\n') for line in f.readlines()]
    except Exception as e:
        print(f"错误：读取Merge数据文件失败：{str(e)}")
        return []
    
    key_mapping = {
        "Mapped Reads Assigned To Intronic Regions Fraction": "Mapped Reads Assigned To Intronic Regions Fraction",
        "Mapped Reads Assigned Antisense To Gene Fraction": "Mapped Reads Assigned Antisense To Gene Fraction"
    }
    
    stats = []
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        parts = stripped.split('\t')
        if len(parts) != 2:
            continue
        stat_name, stat_value = parts[0].strip(), parts[1].strip()
        stat_name = key_mapping.get(stat_name, stat_name)
        processed_value = _process_value(stat_value)
        stats.append((stat_name, processed_value))
        if stat_name in PIE_METRICS and isinstance(processed_value, float):
            merge_pie_data[stat_name] = processed_value
        if stat_name in CARD_METRICS:
            merge_card_data[stat_name] = processed_value
    return stats

def _process_value(value):
    if value.endswith('%'):
        try:
            num = float(value.rstrip('%'))
            return f"{num:.2f}%"
        except ValueError:
            return value
    try:
        float_val = float(value)
        if float_val.is_integer():
            return int(float_val)
        else:
            return float_val
    except ValueError:
        return value

merge_stats = parse_merge_summary(merge_stats_file)
if merge_stats:
    print(f"✓ 解析Merge数据（{len(merge_stats)} 项）")
else:
    print("警告：Merge数据解析失败")

# ---------------------------
# 3.4 解析Merge Barcode数据（08-Sample_Summary/2-Cell_Summary_withSampleID_AllCB.xls）
# 修正：抽样逻辑（保留原有逻辑）
# ---------------------------
merge_barcode_file = "08-Sample_Summary/2-Cell_Summary_withSampleID_AllCB.xls"

def parse_merge_barcode_data(file_path):
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
        if not lines:
            print("警告：Merge Barcode文件为空")
            return []
        header = lines[0].strip().split('\t')
        if 'nUMIunique' not in header:
            print("错误：Merge Barcode文件缺少nUMIunique列")
            return []
        umi_index = header.index('nUMIunique')
        umi_counts = []
        for line in lines[1:]:
            parts = line.strip().split('\t')
            if len(parts) <= umi_index:
                continue
            try:
                umi = int(parts[umi_index])
                umi_counts.append(umi)
            except ValueError:
                print(f"警告：Merge Barcode行 '{line.strip()}' 第{umi_index+1}列非整数，跳过")
                continue
        if not umi_counts:
            print("警告：Merge Barcode无有效数据")
            return []
        
        umi_counts.sort(reverse=True)
        sampled_umi_counts = umi_counts[::SAMPLING_INTERVAL]
        return [{"x": i*SAMPLING_INTERVAL+1, "y": val} for i, val in enumerate(sampled_umi_counts)]
    except FileNotFoundError:
        print(f"错误：{file_path} 不存在")
        return []
    except Exception as e:
        print(f"错误：解析Merge Barcode失败：{str(e)}")
        return []

merge_barcode_rank_data = parse_merge_barcode_data(merge_barcode_file)
if merge_barcode_rank_data:
    print(f"✓ 解析Merge Barcode数据（原始{len(umi_counts)}点，抽样后{len(merge_barcode_rank_data)}点）")

# ---------------------------
# 3.5 解析Sample统计数据（08-Sample_Summary/5-Sample_Summary.xls，Tab分隔）
# 修正：Pandas解析（保留原有逻辑）
# ---------------------------
sample_summary_file = "08-Sample_Summary/5-Sample_Summary.xls"

def parse_sample_summary(file_path):
    if not os.path.exists(file_path):
        print(f"错误：Sample Summary文件 {file_path} 不存在！")
        return {}, []
    try:
        df = pd.read_csv(file_path, sep='\t', header=0, index_col=0)
        df = df.T 
        params = df.index.tolist()
        samples = df.columns.tolist()
        sample_data = {}
        for sample in samples:
            sample_stats = []
            for param in params:
                value = df.loc[param, sample]
                processed_value = _process_value(str(value))
                sample_stats.append((param, processed_value))
            sample_data[sample] = sample_stats
        return sample_data, samples
    except Exception as e:
        print(f"错误：解析Sample Summary失败：{str(e)}")
        return {}, []

def _process_value(value):
    if value.endswith('%'):
        try:
            num = float(value.rstrip('%'))
            return f"{num:.2f}%"
        except ValueError:
            return value
    try:
        float_val = float(value)
        if float_val.is_integer():
            return int(float_val)
        else:
            return float_val
    except ValueError:
        return value

sample_summary_data, sample_samples = parse_sample_summary(sample_summary_file)
if sample_summary_data:
    print(f"✓ 解析Sample数据，发现 {len(sample_samples)} 个样本")
else:
    print("警告：Sample数据解析失败")

# ---------------------------
# 3.6 解析Sample Barcode数据（08-Sample_Summary/2-Cell_Summary_withSampleID_AllCB.xls）
# 修正：抽样逻辑（保留原有逻辑）
# ---------------------------
sample_barcode_file = "08-Sample_Summary/2-Cell_Summary_withSampleID_AllCB.xls"

def parse_sample_barcode_data(file_path):
    if not os.path.exists(file_path):
        print(f"错误：Sample Barcode文件 {file_path} 不存在！")
        return {}
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
        header = lines[0].strip().split('\t')
        sample_id_col = "SampleID"
        if sample_id_col not in header:
            print(f"错误：Sample Barcode文件无 {sample_id_col} 列")
            return {}
        sample_id_index = header.index(sample_id_col)
        umi_index = 16
        sample_umi = {}
        for line in lines[1:]:
            parts = line.strip().split('\t')
            if len(parts) <= max(sample_id_index, umi_index):
                continue
            sample_name = parts[sample_id_index]
            try:
                umi = int(parts[umi_index])
                if sample_name not in sample_umi:
                    sample_umi[sample_name] = []
                sample_umi[sample_name].append(umi)
            except ValueError:
                print(f"警告：Sample Barcode行 '{line.strip()}' 第{umi_index+1}列非整数，跳过")
                continue
        sample_rank = {}
        for sample, umis in sample_umi.items():
            if not umis:
                sample_rank[sample] = []
                continue
            umis.sort(reverse=True)
            sampled_umis = umis[::SAMPLING_INTERVAL]
            sample_rank[sample] = [{"x": i*SAMPLING_INTERVAL+1, "y": val} for i, val in enumerate(sampled_umis)]
        return sample_rank
    except Exception as e:
        print(f"错误：解析Sample Barcode失败：{str(e)}")
        return {}

sample_barcode_rank_data = parse_sample_barcode_data(sample_barcode_file)
if sample_barcode_rank_data:
    print(f"✓ 解析Sample Barcode数据，覆盖 {len(sample_barcode_rank_data)} 个样本")

sample_cards = {
    sample: {stat: val for stat, val in data if stat in CARD_METRICS} 
    for sample, data in sample_summary_data.items()
}

# ---------------------------
# 4. HTML模板（核心优化：修复PDF导出功能）
# ---------------------------
html_template = """
<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <title>10K超高通量单细胞RNA测序分析报告</title>
    <style>
        /* 全局样式（保留基础，优化注释说明模块和按钮样式） */
        body { margin: 0; padding: 0; font-family: Arial, sans-serif; }
        .main-container { max-width: 1200px; margin: 0 auto; padding: 0 50px; }
        header { background: #6A2C70; color: #fff; padding: 12px 0; }
        .header-container { display: flex; align-items: center; gap: 20px; }
        .header-logo { height: 60px; vertical-align: middle; }
        .logo-text { font-size: 24px; font-weight: bold; margin: 0; }
        .subtitle-text { font-size: 16px; margin: 0; }
        .tab-container { display: flex; gap: 10px; margin: 20px 0; }
        .tab-button { 
            padding: 8px 20px; 
            font-size: 16px; 
            background: #eee; 
            border: 1px solid #ccc; 
            border-radius: 4px; 
            cursor: pointer; 
        }
        .tab-button.active { 
            background: #6A2C70; 
            color: #fff; 
            border: none; 
        }
        .tab-content { display: none; }
        .tab-content.active { display: block; }
        .selector { display: flex; align-items: center; margin: 10px 0; }
        .selector span { font-size: 18px; margin-right: 8px; font-weight: bold; }
        select { padding: 8px 15px; font-size: 18px; }
        .content-row { 
            display: flex; 
            align-items: flex-start; 
            gap: 40px; 
            margin: 20px 0; 
        }
        .stats-table, .left-panel { width: 50%; }
        .left-panel { 
            display: flex; 
            flex-direction: column; 
            gap: 20px; 
        }
        .cards-container { 
            display: flex; 
            gap: 25px; 
            flex-wrap: wrap; 
        }
        .stats-table { 
            border-collapse: collapse; 
            font-size: 14px; 
        }
        .stats-table th, .stats-table td { 
            border: 1px solid #CCC; 
            padding: 8px 12px; 
        }
        .stats-table th:first-child, 
        .stats-table td:first-child { 
            text-align: left; 
            white-space: nowrap; 
        }
        .stats-table th:nth-child(2), 
        .stats-table td:nth-child(2) { 
            text-align: right; 
        }
        .stats-table th { 
            background: #F2F2F2; 
            font-weight: bold; 
        }
        .card { 
            background: #fff; 
            border: 1px solid #ccc; 
            border-radius: 8px; 
            padding: 25px; 
            text-align: center; 
            min-width: 200px; 
            box-shadow: 0 4px 8px rgba(0,0,0,0.1); 
            transition: transform 0.2s; 
        }
        .card:hover { transform: scale(1.05); }
        .card-value { 
            font-size: 32px; 
            font-weight: bold; 
            color: #6A2C70; 
            margin-bottom: 12px; 
        }
        .card-label { 
            font-size: 16px; 
            color: #666; 
        }
        /* ===== 注释说明模块优化 ===== */
        .header-with-info { 
            display: flex; 
            align-items: center; 
            justify-content: space-between; 
            padding-right: 40px; 
        }
        .info-details { 
            position: relative; 
            display: inline-flex; /* 改为flex布局，优化图标与文字对齐 */
            align-items: center;
            vertical-align: middle; 
            margin-left: 10px; 
        }
        summary.info-icon { 
            width: 24px; height: 24px; /* 增大图标尺寸，提升可见性 */
            background: #6A2C70; 
            border-radius: 50%; 
            text-align: center; 
            line-height: 24px; /* 调整行高，使问号居中 */
            color: #fff; 
            font-size: 14px; /* 增大字体，更清晰 */
            cursor: pointer; 
            font-weight: bold; 
            border: none; 
            padding: 0; 
            transition: background 0.3s; /* 添加过渡效果，hover时平滑变色 */
        }
        /* 鼠标悬停时提前改变背景色，增强交互反馈 */
        summary.info-icon:hover {
            background: #4A1A4E; 
        }
        /* 展开时旋转图标，增加动态效果 */
        details[open] summary.info-icon { 
            background: #4A1A4E; 
            transform: rotate(45deg); 
            transition: transform 0.3s, background 0.3s; 
        }
        .info-content { 
            width: 600px; /* 调整宽度，更适中 */
            max-height: 400px; /* 增加最大高度，显示更多内容 */
            overflow-y: auto; 
            padding: 16px; /* 增加内边距，提升阅读体验 */
            background: #fff; /* 背景改为白色，更清晰 */
            border: 1px solid #ddd; /* 边框颜色变浅，更柔和 */
            border-radius: 8px; /* 增加圆角，优化外观 */
            font-size: 14px; 
            line-height: 1.6; 
            box-shadow: 0 6px 15px rgba(0,0,0,0.1); /* 增强阴影，提升层次感 */
            position: absolute; 
            left: 50%; /* 水平居中 */
            transform: translateX(-50%); /* 水平居中 */
            top: calc(100% + 10px); /* 距离图标下方10px */
            z-index: 9999; /* 提高层级，避免被覆盖 */
        }
        /* 新增伪元素，添加箭头指示，明确关联图标与说明框 */
        .info-content::before {
            content: '';
            position: absolute;
            top: -10px; /* 箭头位置在内容上方 */
            left: 50%;
            transform: translateX(-50%);
            border-width: 0 10px 10px 10px;
            border-style: solid;
            border-color: transparent transparent #ddd transparent; /* 与边框颜色一致 */
        }
        .info-content::after {
            content: '';
            position: absolute;
            top: -8px; /* 覆盖在伪元素上，形成白色箭头 */
            left: 50%;
            transform: translateX(-50%);
            border-width: 0 10px 10px 10px;
            border-style: solid;
            border-color: transparent transparent #fff transparent; /* 与背景颜色一致 */
        }
        
        /* ===== 图表操作按钮样式 ===== */
        .chart-actions {
            position: absolute;
            top: 10px;
            right: 10px;
            display: flex;
            gap: 5px;
            z-index: 10;
        }
        .action-button { 
            width: 30px;
            height: 30px;
            border-radius: 50%;
            background-color: rgba(255, 255, 255, 0.8);
            display: flex;
            align-items: center;
            justify-content: center;
            cursor: pointer;
            border: 1px solid #ddd;
            transition: all 0.2s;
        }
        .action-button:hover {
            background-color: white;
            box-shadow: 0 2px 5px rgba(0, 0, 0, 0.2);
            transform: translateY(-2px);
        }
        .action-button svg {
            width: 16px;
            height: 16px;
            fill: #6A2C70;
        }
        
        /* ===== 其他样式保留 ===== */
        .no-data { 
            text-align: center; 
            color: #999; 
            font-style: italic; 
            padding: 20px; 
        }
        .svg-chart {
            position: relative;
            width: 100%;
            height: 450px;
            border: 1px solid #eee;
            background-color: white;
            margin-bottom: 20px;
        }
        .axis text { font-family: Arial, sans-serif; font-size: 12px; }
        .axis path, .axis line { fill: none; stroke: #000; shape-rendering: crispEdges; }
        /* 重点修改：移除填充区域样式 */
        .barcode-line { fill: none; stroke: #000; stroke-width: 2px; } /* 曲线改为黑色 */
        .tooltip {
            position: absolute;
            background-color: white;
            border: 1px solid #ccc;
            padding: 8px;
            font-family: Arial, sans-serif;
            font-size: 14px;
            pointer-events: none;
            opacity: 0;
            transition: opacity 0.2s;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }
        .pie-sector { cursor: pointer; transition: transform 0.2s; }
        .pie-sector:hover { transform: scale(1.02); }
        .legend-item { cursor: pointer; }
        .legend-item:hover text { font-weight: bold; }
    </style>
</head>
<body>
    <div class="main-container">
        <div id="global-error" style="display: none;">
            <p>⚠️ 未找到有效数据！</p>
        </div>

        <header>
            <div class="header-container">
                <img src="data:image/png;base64,{{ logo_base64 }}" 
                     alt="万乘基因Logo" class="header-logo">
                <div class="header-text">
                    <div class="logo-text">CellCosmo_UHT Software</div>
                    <div class="subtitle-text">Ultra-High-Throughput Single-Cell RNA Analysis Report</div>
                </div>
            </div>
        </header>

        <div class="tab-container">
            <button class="tab-button active" data-tab="library">Library</button>
            <button class="tab-button" data-tab="merge">MergeLibrary</button>
            <button class="tab-button" data-tab="sample">Sample</button>
        </div>

        <!-- Library选项卡 -->
        <div id="library-content" class="tab-content active">
            <div class="selector">
                <span>Library:</span>  
                <select id="sample-select">
                    {% for sample in library_samples %}
                    <option value="{{ sample }}">{{ sample }}</option>
                    {% endfor %}
                </select>
            </div>
            <div class="content-row">
                <table class="stats-table" id="library-stats-table">
                    <thead>
                        <tr>
                            <th class="header-with-info">
                                <span>Parameter</span>
                                <details class="info-details">
                                    <summary class="info-icon">?</summary>
                                    <div class="info-content">
                                        {% for stat, desc in stat_descriptions.items() %}
                                        <p><strong>{{ stat }}:</strong> {{ desc }}</p>
                                        {% endfor %}
                                    </div>
                                </details>
                            </th>
                            <th>MetricValues</th>
                        </tr>
                    </thead>
                    <tbody></tbody>
                </table>
                <div class="left-panel">
                    <div class="cards-container" id="library-cards-container"></div>
                    <div class="chart-container">
                        <div id="library-barcode-rank-chart" class="svg-chart"></div>
                        <div id="library-region-pie-chart" class="svg-chart"></div>
                    </div>
                </div>
            </div>
            <div id="library-no-data" class="no-data" style="display: none;">
                没有找到该样本的Library数据
            </div>
        </div>

        <!-- Sample选项卡 -->
        <div id="sample-content" class="tab-content">
            <div class="selector">
                <span>Sample:</span>  
                <select id="sample-sample-select">
                    {% for sample in sample_samples %}
                    <option value="{{ sample }}">{{ sample }}</option>
                    {% endfor %}
                </select>
            </div>
            <div class="content-row">
                <table class="stats-table" id="sample-stats-table">
                    <thead>
                        <tr>
                            <th class="header-with-info">
                                <span>Parameter</span>
                                <details class="info-details">
                                    <summary class="info-icon">?</summary>
                                    <div class="info-content">
                                        {% for stat, desc in stat_descriptions.items() %}
                                        <p><strong>{{ stat }}:</strong> {{ desc }}</p>
                                        {% endfor %}
                                    </div>
                                </details>
                            </th>
                            <th>MetricValues</th>
                        </tr>
                    </thead>
                    <tbody></tbody>
                </table>
                <div class="left-panel">
                    <div class="cards-container" id="sample-cards-container"></div>
                    <div class="chart-container">
                        <div id="sample-barcode-rank-chart" class="svg-chart"></div>
                        <div id="sample-region-pie-chart" class="svg-chart"></div>
                    </div>
                </div>
            </div>
            <div id="sample-no-data" class="no-data" style="display: none;">
                没有找到该样本的Sample数据
            </div>
        </div>

        <!-- MergeLibrary选项卡 -->
        <div id="merge-content" class="tab-content">
            <div class="content-row">
                <table class="stats-table" id="merge-stats-table">
                    <thead>
                        <tr>
                            <th class="header-with-info">
                                <span>Parameter</span>
                                <details class="info-details">
                                    <summary class="info-icon">?</summary>
                                    <div class="info-content">
                                        {% for stat, desc in stat_descriptions.items() %}
                                        <p><strong>{{ stat }}:</strong> {{ desc }}</p>
                                        {% endfor %}
                                    </div>
                                </details>
                            </th>
                            <th>MetricValues</th>
                        </tr>
                    </thead>
                    <tbody></tbody>
                </table>
                <div class="left-panel">
                    <div class="cards-container" id="merge-cards-container"></div>
                    <div class="chart-container">
                        <div id="merge-barcode-rank-chart" class="svg-chart"></div>
                        <div id="merge-region-pie-chart" class="svg-chart"></div>
                    </div>
                </div>
            </div>
            <div id="merge-no-data" class="no-data" style="display: none;">
                没有找到MergeLibrary的数据
            </div>
        </div>
    </div>

    <script>
        // 数据传递（保留原有逻辑）
        var libraryData = {{ data_dict|tojson }};
        var libraryBarcode = {{ barcode_rank_data|tojson }};
        var librarySamples = {{ library_samples|tojson }};

        var mergeStats = {{ merge_stats|tojson }};
        var mergeBarcode = {{ merge_barcode_rank_data|tojson }};
        var mergeCard = {{ merge_card_data|tojson }};
        var mergePie = {{ merge_pie_data|tojson }};
        
        var sampleData = {{ sample_summary_data|tojson }};
        var sampleBarcode = {{ sample_barcode_rank_data|tojson }};
        var sampleSamples = {{ sample_samples|tojson }};
        var sampleCards = {{ sample_cards|tojson }};
        
        var cardMetrics = {{ CARD_METRICS|tojson }};
        var pieMetrics = {{ PIE_METRICS|tojson }};
        var statDescriptions = {{ stat_descriptions|tojson }};
        
        // DOM元素（保留原有逻辑）
        var tabButtons = document.querySelectorAll('.tab-button');
        var tabContents = document.querySelectorAll('.tab-content');
        var librarySelect = document.getElementById('sample-select');
        var sampleSelect = document.getElementById('sample-sample-select');
        
        var libraryTable = document.getElementById('library-stats-table').getElementsByTagName('tbody')[0];
        var libraryCards = document.getElementById('library-cards-container');
        var libraryBarcodeChart = document.getElementById('library-barcode-rank-chart');
        var libraryPieChart = document.getElementById('library-region-pie-chart');
        var libraryNoData = document.getElementById('library-no-data');
        
        var mergeTable = document.getElementById('merge-stats-table').getElementsByTagName('tbody')[0];
        var mergeCards = document.getElementById('merge-cards-container');
        var mergeBarcodeChart = document.getElementById('merge-barcode-rank-chart');
        var mergePieChart = document.getElementById('merge-region-pie-chart');
        var mergeNoData = document.getElementById('merge-no-data');

        var sampleTable = document.getElementById('sample-stats-table').getElementsByTagName('tbody')[0];
        var sampleCardsContainer = document.getElementById('sample-cards-container');
        var sampleBarcodeChart = document.getElementById('sample-barcode-rank-chart');
        var samplePieChart = document.getElementById('sample-region-pie-chart');
        var sampleNoData = document.getElementById('sample-no-data');
        
        // 工具函数（保留原有逻辑）
        function toRadians(degrees) {
            return degrees * Math.PI / 180;
        }
        function createSvgElement(tag, attrs = {}) {
            const elem = document.createElementNS("http://www.w3.org/2000/svg", tag);
            Object.entries(attrs).forEach(([k, v]) => elem.setAttribute(k, v));
            return elem;
        }
        function calculateLogTicks(min, max, count = 5) {
            const minExp = Math.floor(Math.log10(min));
            const maxExp = Math.ceil(Math.log10(max));
            const ticks = [];
            for (let exp = minExp; exp <= maxExp; exp++) {
                for (let i = 1; i <= 9; i++) {
                    const tick = i * Math.pow(10, exp);
                    if (tick >= min && tick <= max) {
                        ticks.push(tick);
                        if (ticks.length >= count) break;
                    }
                }
                if (ticks.length >= count) break;
            }
            return ticks;
        }
        
        // 新增：将SVG转换为PDF的函数（使用纯JavaScript实现）
        function svgToPdf(svgElement, fileName) {
            return new Promise((resolve, reject) => {
                try {
                    // 创建临时canvas
                    const canvas = document.createElement('canvas');
                    const ctx = canvas.getContext('2d');
                    
                    // 获取SVG的尺寸
                    const svgRect = svgElement.getBoundingClientRect();
                    const width = svgRect.width;
                    const height = svgRect.height;
                    
                    // 设置canvas尺寸为SVG尺寸的2倍（提高分辨率）
                    const scale = 2; // 缩放因子
                    canvas.width = width * scale;
                    canvas.height = height * scale;
                    
                    // 将SVG转换为图像
                    const svgString = new XMLSerializer().serializeToString(svgElement);
                    const svgBlob = new Blob([svgString], { type: 'image/svg+xml;charset=utf-8' });
                    const svgUrl = URL.createObjectURL(svgBlob);
                    
                    const img = new Image();
                    img.onload = function() {
                        // 绘制图像到canvas（缩放以提高分辨率）
                        ctx.scale(scale, scale);
                        ctx.fillStyle = 'white';
                        ctx.fillRect(0, 0, width, height);
                        ctx.drawImage(img, 0, 0, width, height);
                        
                        // 从canvas获取数据URL
                        const imgData = canvas.toDataURL('image/png');
                        
                        // 创建PDF文档
                        const pdf = new jsPDF('landscape', 'pt', [width, height]);
                        
                        // 添加图像到PDF
                        pdf.addImage(imgData, 'PNG', 0, 0, width, height);
                        
                        // 保存PDF
                        pdf.save(fileName);
                        
                        // 清理资源
                        URL.revokeObjectURL(svgUrl);
                        
                        resolve();
                    };
                    
                    img.onerror = function() {
                        reject(new Error('Failed to load SVG image'));
                        URL.revokeObjectURL(svgUrl);
                    };
                    
                    img.src = svgUrl;
                } catch (error) {
                    reject(error);
                }
            });
        }
        
        // 动态加载jsPDF库
        function loadJsPdf() {
            return new Promise((resolve, reject) => {
                if (window.jsPDF) {
                    resolve();
                    return;
                }
                
                const script = document.createElement('script');
                script.src = 'https://cdnjs.cloudflare.com/ajax/libs/jspdf/2.5.1/jspdf.umd.min.js';
                script.onload = () => {
                    if (window.jspdf && window.jspdf.jsPDF) {
                        window.jsPDF = window.jspdf.jsPDF;
                        resolve();
                    } else {
                        reject(new Error('jsPDF library failed to load'));
                    }
                };
                script.onerror = () => reject(new Error('Failed to load jsPDF library'));
                document.head.appendChild(script);
            });
        }
        
        // Barcode Rank图渲染（核心修改：移除填充区域，曲线改为黑色）
        function renderBarcodeChart(sampleName, container, dataPoints, isMerge = false) {
            if (dataPoints.length === 0) {
                container.innerHTML = '<div class="no-data">无Barcode数据</div>';
                return;
            }
            container.innerHTML = '';
            const width = container.clientWidth;
            const height = 450;
            const margin = { top: 40, right: 20, bottom: 60, left: 70 };
            const innerWidth = width - margin.left - margin.right;
            const innerHeight = height - margin.top - margin.bottom;
            const svg = createSvgElement('svg', {
                class: 'svg-chart-svg',
                viewBox: `0 0 ${width} ${height}`
            });
            container.appendChild(svg);
            
            // 添加图表操作按钮
            const actions = document.createElement('div');
            actions.className = 'chart-actions';
            actions.innerHTML = `
                <div class="action-button" id="capture-pdf" title="下载PDF矢量图">
                    <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24">
                        <path d="M19 9h-4V3H9v6H5l7 7 7-7zM5 18v2h14v-2H5z"/>
                    </svg>
                </div>
                <div class="action-button" id="home-btn" title="重置视图">
                    <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24">
                        <path d="M10 20v-6h4v6h5v-8h3L12 3 2 12h3v8h5z"/>
                    </svg>
                </div>
            `;
            container.appendChild(actions);
            
            const g = createSvgElement('g', {
                transform: `translate(${margin.left}, ${margin.top})`
            });
            svg.appendChild(g);
            const xValues = dataPoints.map(p => p.x);
            const yValues = dataPoints.map(p => p.y);
            const xMin = Math.max(1, Math.min(...xValues));
            const xMax = Math.max(...xValues);
            const yMin = Math.max(1, Math.min(...yValues));
            const yMax = Math.max(...yValues);
            function logScale(val, min, max, rMin, rMax) {
                const logVal = Math.log(val);
                const logMin = Math.log(min);
                const logMax = Math.log(max);
                return rMin + (logVal - logMin) * (rMax - rMin) / (logMax - logMin);
            }
            const xTicks = calculateLogTicks(xMin, xMax, 5);
            const yTicks = calculateLogTicks(yMin, yMax, 5);
            yTicks.forEach(tick => {
                const y = innerHeight - logScale(tick, yMin, yMax, 0, innerHeight);
                const line = createSvgElement('line', {
                    x1: 0, y1: y, x2: innerWidth, y2: y,
                    stroke: '#f0f0f0', // 改为浅灰色
                    'stroke-width': 1
                });
                g.appendChild(line);
            });
            xTicks.forEach(tick => {
                const x = logScale(tick, xMin, xMax, 0, innerWidth);
                const line = createSvgElement('line', {
                    x1: x, y1: 0, x2: x, y2: innerHeight,
                    stroke: '#f0f0f0', // 改为浅灰色
                    'stroke-width': 1
                });
                g.appendChild(line);
            });
            const lineData = dataPoints.map(p => ({
                x: logScale(p.x, xMin, xMax, 0, innerWidth),
                y: innerHeight - logScale(p.y, yMin, yMax, 0, innerHeight)
            }));
            
            // 核心修改：只绘制曲线，不绘制填充区域
            const linePath = lineData.reduce((p, pt, i) => `${p} ${i===0?'M':'L'} ${pt.x} ${pt.y}`, '');
            const line = createSvgElement('path', { 
                d: linePath, 
                fill: "none",  // 确保无填充
                stroke: "#000", // 曲线改为黑色
                'stroke-width': 2 
            });
            g.appendChild(line);
            
            const xAxis = createSvgElement('g', {
                transform: `translate(0, ${innerHeight})`,
                class: 'axis'
            });
            g.appendChild(xAxis);
            xTicks.forEach(tick => {
                const x = logScale(tick, xMin, xMax, 0, innerWidth);
                const tickGroup = createSvgElement('g');
                const tickLine = createSvgElement('line', {
                    x1: x, y1: 0, x2: x, y2: 5,
                    stroke: '#000'
                });
                tickGroup.appendChild(tickLine);
                const tickText = createSvgElement('text', {
                    x: x, y: 20,
                    'text-anchor': 'middle',
                    fill: '#333'
                });
                tickText.textContent = tick >= 1000 ? `${tick/1000}k` : tick;
                tickGroup.appendChild(tickText);
                xAxis.appendChild(tickGroup);
            });
            const xTitle = createSvgElement('text', {
                x: innerWidth/2, y: 40,
                'text-anchor': 'middle',
                'font-weight': 'bold',
                'font-size': '14px'
            });
            xTitle.textContent = 'Barcode Rank (log₁₀)';
            xAxis.appendChild(xTitle);
            const yAxis = createSvgElement('g', { class: 'axis' });
            g.appendChild(yAxis);
            yTicks.forEach(tick => {
                const y = innerHeight - logScale(tick, yMin, yMax, 0, innerHeight);
                const tickGroup = createSvgElement('g');
                const tickLine = createSvgElement('line', {
                    x1: 0, y1: y, x2: -5, y2: y,
                    stroke: '#000'
                });
                tickGroup.appendChild(tickLine);
                const tickText = createSvgElement('text', {
                    x: -10, y: y + 4,
                    'text-anchor': 'end',
                    fill: '#333'
                });
                tickText.textContent = tick >= 1000 ? `${tick/1000}k` : tick;
                tickGroup.appendChild(tickText);
                yAxis.appendChild(tickGroup);
            });
            const yTitle = createSvgElement('text', {
                x: -innerHeight/2, y: -40,
                'text-anchor': 'middle',
                'font-weight': 'bold',
                'font-size': '14px',
                transform: 'rotate(-90)'
            });
            yTitle.textContent = 'UMI Counts (log₁₀)';
            yAxis.appendChild(yTitle);
            const title = createSvgElement('text', {
                x: innerWidth/2, y: -10,
                'text-anchor': 'middle',
                'font-weight': 'bold',
                'font-size': '16px'
            });
            title.textContent = `${sampleName} Barcode Rank`;
            g.appendChild(title);
            const tooltip = document.createElement('div');
            tooltip.className = 'tooltip';
            document.body.appendChild(tooltip);
            const interactive = createSvgElement('rect', {
                x: 0, y: 0,
                width: innerWidth, height: innerHeight,
                fill: 'transparent',
                'pointer-events': 'all'
            });
            g.appendChild(interactive);
            
            // 图表平移和缩放功能
            let isDragging = false;
            let startX, startY;
            let currentTransform = { x: 0, y: 0, scale: 1 };
            let originalTransform = { x: 0, y: 0, scale: 1 };
            
            const transformGroup = g; // 变换应用到g元素上
            
            function applyTransform() {
                transformGroup.setAttribute('transform', 
                    `translate(${margin.left + currentTransform.x}, ${margin.top + currentTransform.y}) scale(${currentTransform.scale})`
                );
            }
            
            interactive.addEventListener('mousedown', (e) => {
                isDragging = true;
                startX = e.clientX;
                startY = e.clientY;
                document.addEventListener('mousemove', onMouseMove);
                document.addEventListener('mouseup', onMouseUp);
            });
            
            function onMouseMove(e) {
                if (!isDragging) return;
                const dx = e.clientX - startX;
                const dy = e.clientY - startY;
                currentTransform.x += dx;
                currentTransform.y += dy;
                applyTransform();
                startX = e.clientX;
                startY = e.clientY;
            }
            
            function onMouseUp() {
                isDragging = false;
                document.removeEventListener('mousemove', onMouseMove);
                document.removeEventListener('mouseup', onMouseUp);
            }
            
            // 鼠标滚轮缩放
            interactive.addEventListener('wheel', (e) => {
                e.preventDefault();
                const scaleFactor = e.deltaY < 0 ? 1.1 : 0.9;
                const mouseX = e.clientX;
                const mouseY = e.clientY;
                
                // 计算鼠标在SVG中的位置
                const rect = svg.getBoundingClientRect();
                const svgX = (mouseX - rect.left - margin.left - currentTransform.x) / currentTransform.scale;
                const svgY = (mouseY - rect.top - margin.top - currentTransform.y) / currentTransform.scale;
                
                // 应用缩放
                currentTransform.scale *= scaleFactor;
                currentTransform.x = mouseX - rect.left - svgX * currentTransform.scale;
                currentTransform.y = mouseY - rect.top - svgY * currentTransform.scale;
                
                applyTransform();
            });
            
            // PDF按钮点击事件
            actions.querySelector('#capture-pdf').addEventListener('click', () => {
                const fileName = `${sampleName || 'barcode'}_rank.pdf`;
                
                // 加载jsPDF库（如果需要）
                loadJsPdf()
                    .then(() => {
                        return svgToPdf(svg, fileName);
                    })
                    .then(() => {
                        showNotification(`已保存 ${fileName}`);
                    })
                    .catch(err => {
                        showNotification(`保存失败: ${err.message}`);
                    });
            });
            
            // Home按钮点击事件
            actions.querySelector('#home-btn').addEventListener('click', () => {
                currentTransform = { ...originalTransform };
                applyTransform();
                showNotification('已恢复原始视图');
            });
            
            interactive.addEventListener('mousemove', (e) => {
                const rect = svg.getBoundingClientRect();
                const mouseX = e.clientX - rect.left - margin.left - currentTransform.x;
                let closest = null;
                let minDist = Infinity;
                lineData.forEach((pt, i) => {
                    const dist = Math.abs(pt.x - mouseX / currentTransform.scale);
                    if (dist < minDist) {
                        minDist = dist;
                        closest = dataPoints[i];
                    }
                });
                if (closest) {
                    tooltip.style.opacity = '1';
                    tooltip.style.left = `${e.clientX + 10}px`;
                    tooltip.style.top = `${e.clientY + 10}px`;
                    tooltip.innerHTML = `
                        <div><strong>Rank:</strong> ${closest.x.toLocaleString()}</div>
                        <div><strong>UMI:</strong> ${closest.y.toLocaleString()}</div>
                        <div><strong>log₁₀(Rank):</strong> ${Math.log10(closest.x).toFixed(2)}</div>
                        <div><strong>log₁₀(UMI):</strong> ${Math.log10(closest.y).toFixed(2)}</div>
                    `;
                }
            });
            interactive.addEventListener('mouseout', () => {
                tooltip.style.opacity = '0';
            });
        }
        
        // 饼图渲染（保留原有逻辑，增加PDF导出）
        function renderPieChart(exonic, intronic, intergenic, antisense, container, title) {
            container.innerHTML = '';
            const regions = [
                { name: 'Exonic Regions', value: exonic, color: '#FF6384' },
                { name: 'Intronic Regions', value: intronic, color: '#36A2EB' },
                { name: 'Intergenic Regions', value: intergenic, color: '#CC65FE' },
                { name: 'Antisense To Gene', value: antisense, color: '#FFCE56' }
            ].filter(r => r.value > 0);
            if (regions.length === 0) {
                container.innerHTML = '<div class="no-data">无区域分布数据</div>';
                return;
            }
            const width = container.clientWidth;
            const height = 450;
            const margin = { top: 40, right: 200, bottom: 20, left: 20 };
            const innerWidth = width - margin.left - margin.right;
            const innerHeight = height - margin.top - margin.bottom;
            const radius = Math.min(innerWidth, innerHeight) / 2;
            const svg = createSvgElement('svg', {
                class: 'svg-chart-svg',
                viewBox: `0 0 ${width} ${height}`
            });
            container.appendChild(svg);
            
            // 添加图表操作按钮
            const actions = document.createElement('div');
            actions.className = 'chart-actions';
            actions.innerHTML = `
                <div class="action-button" id="capture-pdf" title="下载PDF矢量图">
                    <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24">
                        <path d="M19 9h-4V3H9v6H5l7 7 7-7zM5 18v2h14v-2H5z"/>
                    </svg>
                </div>
                <div class="action-button" id="home-btn" title="重置视图">
                    <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24">
                        <path d="M10 20v-6h4v6h5v-8h3L12 3 2 12h3v8h5z"/>
                    </svg>
                </div>
            `;
            container.appendChild(actions);
            
            const g = createSvgElement('g', {
                transform: `translate(${margin.left + innerWidth/2}, ${margin.top + innerHeight/2})`
            });
            svg.appendChild(g);
            const total = regions.reduce((s, r) => s + r.value, 0);
            let startAngle = 0;
            regions.forEach((r, i) => {
                const sliceAngle = (r.value / total) * 360;
                const endAngle = startAngle + sliceAngle;
                const startX = Math.cos(toRadians(startAngle)) * radius;
                const startY = Math.sin(toRadians(startAngle)) * radius;
                const endX = Math.cos(toRadians(endAngle)) * radius;
                const endY = Math.sin(toRadians(endAngle)) * radius;
                const largeArc = sliceAngle > 180 ? 1 : 0;
                const pathData = [
                    `M 0 0`,
                    `L ${startX} ${startY}`,
                    `A ${radius} ${radius} 0 ${largeArc} 1 ${endX} ${endY}`,
                    `Z`
                ].join(' ');
                const sector = createSvgElement('path', {
                    d: pathData,
                    fill: r.color,
                    stroke: 'white',
                    'stroke-width': 2,
                    class: 'pie-sector',
                    'data-value': r.value,
                    'data-percentage': ((r.value/total)*100).toFixed(2),
                    'data-name': r.name
                });
                g.appendChild(sector);
                if (sliceAngle > 10) {
                    const centerAngle = startAngle + sliceAngle/2;
                    const centerX = Math.cos(toRadians(centerAngle)) * (radius * 0.7);
                    const centerY = Math.sin(toRadians(centerAngle)) * (radius * 0.7);
                    const label = createSvgElement('text', {
                        x: centerX, y: centerY,
                        'text-anchor': 'middle',
                        'dominant-baseline': 'middle',
                        fill: 'white',
                        'font-weight': 'bold',
                        'font-size': '12px'
                    });
                    label.textContent = `${((r.value/total)*100).toFixed(1)}%`;
                    g.appendChild(label);
                }
                startAngle = endAngle;
            });
            const legend = createSvgElement('g', {
                transform: `translate(${margin.left + innerWidth + 20}, ${margin.top})`
            });
            svg.appendChild(legend);
            regions.forEach((r, i) => {
                const item = createSvgElement('g', {
                    transform: `translate(0, ${i*25})`,
                    class: 'legend-item'
                });
                legend.appendChild(item);
                const box = createSvgElement('rect', {
                    x: 0, y: 0,
                    width: 15, height: 15,
                    fill: r.color,
                    stroke: 'white',
                    'stroke-width': 1
                });
                item.appendChild(box);
                const text = createSvgElement('text', {
                    x: 25, y: 12,
                    'font-size': '12px'
                });
                text.textContent = `${r.name}: ${((r.value/total)*100).toFixed(1)}%`;
                item.appendChild(text);
                item.addEventListener('click', () => {
                    const sector = g.querySelector(`.pie-sector[data-name="${r.name}"]`);
                    if (sector.getAttribute('opacity') === '0') {
                        sector.setAttribute('opacity', '1');
                        text.setAttribute('text-decoration', 'none');
                    } else {
                        sector.setAttribute('opacity', '0');
                        text.setAttribute('text-decoration', 'line-through');
                    }
                });
            });
            const titleElem = createSvgElement('text', {
                x: 40, y: -innerHeight/2 - 10,
                'text-anchor': 'middle',
                'font-weight': 'bold',
                'font-size': '16px'
            });
            titleElem.textContent = title || 'Mapped Reads Region Distribution';
            g.appendChild(titleElem);
            const tooltip = document.createElement('div');
            tooltip.className = 'tooltip';
            document.body.appendChild(tooltip);
            
            // 图表平移和缩放功能
            let isDragging = false;
            let startX, startY;
            let currentTransform = { x: 0, y: 0, scale: 1 };
            let originalTransform = { x: 0, y: 0, scale: 1 };
            
            const transformGroup = g; // 变换应用到g元素上
            
            function applyTransform() {
                transformGroup.setAttribute('transform', 
                    `translate(${margin.left + innerWidth/2 + currentTransform.x}, ${margin.top + innerHeight/2 + currentTransform.y}) scale(${currentTransform.scale})`
                );
            }
            
            svg.addEventListener('mousedown', (e) => {
                if (e.target.classList.contains('pie-sector') || e.target.closest('.legend-item')) return; // 不捕获扇区和图例的点击
                isDragging = true;
                startX = e.clientX;
                startY = e.clientY;
                document.addEventListener('mousemove', onMouseMove);
                document.addEventListener('mouseup', onMouseUp);
            });
            
            function onMouseMove(e) {
                if (!isDragging) return;
                const dx = e.clientX - startX;
                const dy = e.clientY - startY;
                currentTransform.x += dx;
                currentTransform.y += dy;
                applyTransform();
                startX = e.clientX;
                startY = e.clientY;
            }
            
            function onMouseUp() {
                isDragging = false;
                document.removeEventListener('mousemove', onMouseMove);
                document.removeEventListener('mouseup', onMouseUp);
            }
            
            // 鼠标滚轮缩放
            svg.addEventListener('wheel', (e) => {
                if (e.target.classList.contains('pie-sector') || e.target.closest('.legend-item')) return; // 不捕获扇区和图例的滚轮事件
                e.preventDefault();
                const scaleFactor = e.deltaY < 0 ? 1.1 : 0.9;
                const mouseX = e.clientX;
                const mouseY = e.clientY;
                
                // 计算鼠标在SVG中的位置
                const rect = svg.getBoundingClientRect();
                const svgX = (mouseX - rect.left - margin.left - innerWidth/2 - currentTransform.x) / currentTransform.scale;
                const svgY = (mouseY - rect.top - margin.top - innerHeight/2 - currentTransform.y) / currentTransform.scale;
                
                // 应用缩放
                currentTransform.scale *= scaleFactor;
                currentTransform.x = mouseX - rect.left - margin.left - innerWidth/2 - svgX * currentTransform.scale;
                currentTransform.y = mouseY - rect.top - margin.top - innerHeight/2 - svgY * currentTransform.scale;
                
                applyTransform();
            });
            
            // PDF按钮点击事件
            actions.querySelector('#capture-pdf').addEventListener('click', () => {
                const fileName = title ? title.replace(/ /g, '_') + '.pdf' : 'pie_chart.pdf';
                
                // 加载jsPDF库（如果需要）
                loadJsPdf()
                    .then(() => {
                        return svgToPdf(svg, fileName);
                    })
                    .then(() => {
                        showNotification(`已保存 ${fileName}`);
                    })
                    .catch(err => {
                        showNotification(`保存失败: ${err.message}`);
                    });
            });
            
            // Home按钮点击事件
            actions.querySelector('#home-btn').addEventListener('click', () => {
                currentTransform = { ...originalTransform };
                applyTransform();
                showNotification('已恢复原始视图');
            });
            
            g.querySelectorAll('.pie-sector').forEach(sector => {
                sector.addEventListener('mouseover', (e) => {
                    const rect = svg.getBoundingClientRect();
                    const val = sector.getAttribute('data-value');
                    const perc = sector.getAttribute('data-percentage');
                    const name = sector.getAttribute('data-name');
                    tooltip.style.opacity = '1';
                    tooltip.style.left = `${e.clientX + 10}px`;
                    tooltip.style.top = `${e.clientY + 10}px`;
                    tooltip.innerHTML = `
                        <div><strong>${name}</strong></div>
                        <div>值: ${parseFloat(val).toLocaleString()}</div>
                        <div>百分比: ${perc}%</div>
                    `;
                    sector.setAttribute('transform', 'scale(1.03)');
                });
                sector.addEventListener('mouseout', () => {
                    tooltip.style.opacity = '0';
                    sector.setAttribute('transform', 'scale(1)');
                });
            });
        }

        // Library渲染函数（保留原有逻辑）
        function renderLibrary(sampleName) {
            libraryTable.innerHTML = '';
            libraryCards.innerHTML = '';
            libraryBarcodeChart.innerHTML = '';
            libraryPieChart.innerHTML = '';
            libraryNoData.style.display = 'none';
            const stats = libraryData[sampleName];
            if (!stats || stats.length === 0) {
                libraryNoData.style.display = 'block';
                return;
            }
            stats.forEach(([stat, val]) => {
                const row = document.createElement('tr');
                row.innerHTML = `<td>${stat}</td><td>${val}</td>`;
                libraryTable.appendChild(row);
            });
            const cards = {};
            stats.forEach(([stat, val]) => {
                if (cardMetrics.includes(stat)) {
                    cards[stat] = val;
                }
            });
            cardMetrics.forEach(metric => {
                const val = cards[metric] || 'N/A';
                const div = document.createElement('div');
                div.className = 'card';
                div.innerHTML = `
                    <div class="card-value">${val}</div>
                    <div class="card-label">${metric}</div>
                `;
                libraryCards.appendChild(div);
            });
            renderBarcodeChart(sampleName, libraryBarcodeChart, libraryBarcode[sampleName] || []);
            const pieData = { exonic: 0, intronic: 0, intergenic: 0, antisense: 0 };
            stats.forEach(([stat, val]) => {
                let numVal = parseFloat(val.replace('%', '')) || 0;
                if (stat.includes(pieMetrics[0])) pieData.exonic = numVal;
                if (stat.includes(pieMetrics[1])) pieData.intronic = numVal;
                if (stat.includes(pieMetrics[2])) pieData.intergenic = numVal;
                if (stat.includes(pieMetrics[3])) pieData.antisense = numVal;
            });
            renderPieChart(
                pieData.exonic, 
                pieData.intronic, 
                pieData.intergenic, 
                pieData.antisense, 
                libraryPieChart, 
                `${sampleName} Mapped Reads Region Distribution`
            );
        }
        
        // MergeLibrary渲染函数（保留原有逻辑）
        function renderMerge() {
            mergeTable.innerHTML = '';
            mergeCards.innerHTML = '';
            mergeBarcodeChart.innerHTML = '';
            mergePieChart.innerHTML = '';
            mergeNoData.style.display = 'none';
            if (!mergeStats || mergeStats.length === 0) {
                mergeNoData.style.display = 'block';
                return;
            }
            mergeStats.forEach(([stat, val]) => {
                const row = document.createElement('tr');
                row.innerHTML = `<td>${stat}</td><td>${val}</td>`;
                mergeTable.appendChild(row);
            });
            cardMetrics.forEach(metric => {
                const val = mergeCard[metric] || 'N/A';
                const div = document.createElement('div');
                div.className = 'card';
                div.innerHTML = `
                    <div class="card-value">${val}</div>
                    <div class="card-label">${metric}</div>
                `;
                mergeCards.appendChild(div);
            });
            renderBarcodeChart('MergeLibrary', mergeBarcodeChart, mergeBarcode || [], "mergeBarcode");
            renderPieChart(
                mergePie[pieMetrics[0]] || 0,
                mergePie[pieMetrics[1]] || 0,
                mergePie[pieMetrics[2]] || 0,
                mergePie[pieMetrics[3]] || 0,
                mergePieChart,
                'MergeLibrary Mapped Reads Region Distribution',
                "mergePie"
            );
        }

        // Sample渲染函数（保留原有逻辑）
        function renderSample(sampleName) {
            sampleTable.innerHTML = '';
            sampleCardsContainer.innerHTML = '';
            sampleBarcodeChart.innerHTML = '';
            samplePieChart.innerHTML = '';
            sampleNoData.style.display = 'none';
            const stats = sampleData[sampleName];
            if (!stats || stats.length === 0) {
                sampleNoData.style.display = 'block';
                return;
            }
            
            console.log(`[Sample ${sampleName}] 统计项列表:`, stats.map(([stat, _]) => stat));
            stats.forEach(([stat, val]) => {
                const row = document.createElement('tr');
                row.innerHTML = `<td>${stat}</td><td>${val}</td>`;
                sampleTable.appendChild(row);
            });
            const cards = sampleCards[sampleName] || {};
            cardMetrics.forEach(metric => {
                const val = cards[metric] || 'N/A';
                const div = document.createElement('div');
                div.className = 'card';
                div.innerHTML = `
                    <div class="card-value">${val}</div>
                    <div class="card-label">${metric}</div>
                `;
                sampleCardsContainer.appendChild(div);
            });
            renderBarcodeChart(sampleName, sampleBarcodeChart, sampleBarcode[sampleName] || [], "sampleBarcode");

            const pieData = { exonic: 0, intronic: 0, intergenic: 0, antisense: 0 };
            stats.forEach(([stat, val]) => {
                let numVal = typeof val === 'number' ? val : 0;
                if (stat === pieMetrics[0]) pieData.exonic = numVal;
                if (stat === pieMetrics[1]) pieData.intronic = numVal;
                if (stat === pieMetrics[2]) pieData.intergenic = numVal;
                if (stat === pieMetrics[3]) pieData.antisense = numVal;
            });

            renderPieChart(
                pieData.exonic, 
                pieData.intronic, 
                pieData.intergenic, 
                pieData.antisense, 
                samplePieChart, 
                `${sampleName} Mapped Reads Region Distribution`,
                "samplePie"
            );
        }

        
        // 事件监听（保留原有逻辑）
        tabButtons.forEach(button => {
            button.addEventListener('click', () => {
                const tabId = button.getAttribute('data-tab');
                tabButtons.forEach(btn => btn.classList.remove('active'));
                button.classList.add('active');
                tabContents.forEach(content => {
                    content.classList.remove('active');
                    if (content.id === `${tabId}-content`) {
                        content.classList.add('active');
                    }
                });
                if (tabId === 'library') {
                    renderLibrary(librarySelect.value);
                } else if (tabId === 'merge') {
                    renderMerge();
                } else if (tabId === 'sample') {
                    renderSample(sampleSelect.value);
                }
            });
        });
        
        librarySelect.addEventListener('change', (e) => {
            renderLibrary(e.target.value);
        });
        
        sampleSelect.addEventListener('change', (e) => {
            renderSample(e.target.value);
        });
        
        // 初始化渲染（保留原有逻辑）
        if (librarySamples.length > 0) {
            renderLibrary(librarySamples[0]);
        } else {
            document.getElementById('global-error').style.display = 'block';
        }
        
        // 响应式处理（保留原有逻辑）
        window.addEventListener('resize', () => {
            const activeTab = document.querySelector('.tab-button.active').getAttribute('data-tab');
            if (activeTab === 'library') {
                renderLibrary(librarySelect.value);
            } else if (activeTab === 'sample') {
                renderSample(sampleSelect.value);
            } else if (activeTab === 'merge') {
                renderMerge();
            }
        });
    </script>
</body>
</html>
"""

# ---------------------------
# 5. 渲染并保存HTML报告（保留原有逻辑）
# ---------------------------
def generate_report(html_template, output_file="CellCosmo_UHT_Analysis_Report.html"):
    try:
        library_samples = list(data_dict.keys()) if data_dict else []
        template = Template(html_template)
        rendered_html = template.render(
            data_dict=data_dict,
            barcode_rank_data=barcode_rank_data,
            sample_summary_data=sample_summary_data,
            sample_barcode_rank_data=sample_barcode_rank_data,
            sample_samples=sample_samples,
            sample_cards=sample_cards,
            merge_stats=merge_stats,
            merge_card_data=merge_card_data,
            merge_pie_data=merge_pie_data,
            merge_barcode_rank_data=merge_barcode_rank_data,
            CARD_METRICS=CARD_METRICS,
            PIE_METRICS=PIE_METRICS,
            stat_descriptions=stat_descriptions,
            logo_base64=logo_base64,
            library_samples=library_samples
        )
        with open(output_file, 'w') as f:
            f.write(rendered_html)
        print(f"✓ 报告已成功生成：{output_file}")
        return output_file
    except Exception as e:
        print(f"错误：生成报告失败：{str(e)}")
        return None

# ---------------------------
# 6. 主程序入口（保留原有逻辑）
# ---------------------------
if __name__ == "__main__":
    output_file = generate_report(html_template)
    if output_file:
        print("请在浏览器中打开生成的HTML文件查看报告")