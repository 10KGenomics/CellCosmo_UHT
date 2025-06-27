# -*- coding: utf-8 -*-
"""
AddSampleID.py

功能：
    从Cell_Summary.xls（纯文本、Tab分隔）文件中读取细胞级统计数据，
    根据CB列（格式：<seg1>_<seg2>_<seg3>）的第二段barcode，将每个细胞分配到任意数量的Sample组。
    barcode库由外部文件Barcode2.list给出，范围和样本名称通过命令行参数灵活指定。

用法：
    python AddSampleID.py \
        --inputfile 2-Cell_Summary.xls \
        --barcodelist Barcode2.list \
        --splitCB 1-64 65-128 129-192 \
        --splitSample SampleA SampleB SampleC \
        --output 2-Cell_Summary_withSampleID.xls

    上述示例中：
      - 索引1-64映射到SampleA；65-128映射到SampleB；129-192映射到SampleC。
    如果只需要一个样本，可提供单个范围和名称；两个样本则提供两个，以此类推。

原理：
    1. 使用argparse读取并验证可变数量的--splitCB和--splitSample参数数量一致；
    2. 加载Barcode2.list，生成按顺序排列的barcode列表；
    3. 解析每个范围字符串"start-end"为整数区间，并检查不越界；
    4. 根据每个区间和对应样本名称构建barcode->SampleID映射；
    5. 读取Cell_Summary表格，拆分CB列提取第二段barcode，映射SampleID；
    6. 将SampleID列追加到DataFrame，并以Tab分隔格式写出。

作者：
    自动生成脚本，详细注释说明每一步。
"""

import argparse
import sys
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description='Step9.AddSampleID: 根据CB列第二段barcode给每行添加SampleID（可变组数）')
    parser.add_argument(
        '--inputfile', '-i', required=True,
        help='输入的Cell_Summary.xls（Tab分隔）文件路径')
    parser.add_argument(
        '--barcodelist', '-b', required=True,
        help='包含barcode的文件，每行一个barcode，顺序定义样本分组')
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
        help='与splitCB数量对应的样本名称列表，顺序需严格匹配'
    )
    parser.add_argument(
        '--output', '-o', required=True,
        help='输出文件路径，带SampleID列的Cell_Summary_withSampleID.xls')
    return parser.parse_args()


def load_barcodes(path):
    """
    读取barcode列表文件，返回barcode_list（保持顺序）
    """
    with open(path) as f:
        barcodes = [line.strip() for line in f if line.strip()]
    if not barcodes:
        print(f"ERROR: barcode列表文件{path}为空", file=sys.stderr)
        sys.exit(1)
    return barcodes


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
            print(f"ERROR: 无效的等差数列格式 '{spec}': {e}", file=sys.stderr)
            sys.exit(1)
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
                indices.extend(range(a - 1, b))  # range是左闭右开，b不转换因为原格式包含end
            except:
                print(f"ERROR: 无效的范围格式 '{part}'，应为1-based的start-end", file=sys.stderr)
                sys.exit(1)
        else:
            # 单个索引
            try:
                idx = int(part)
                if idx < 1 or idx > total_length:
                    raise ValueError
                indices.append(idx - 1)  # 转换为0-based
            except:
                print(f"ERROR: 无效的单个索引 '{part}'，应为1-based整数", file=sys.stderr)
                sys.exit(1)

    # 去重并保持顺序
    seen = set()
    return [x for x in indices if x not in seen and not seen.add(x)]


def build_mapping(barcodes, specs, samples):
    """
    根据specs和samples构建barcode->SampleID映射字典
    specs: list of 索引规范字符串
    samples: list of样本名称，长度与specs相同
    返回 dict: barcode(str) -> sample_id(str)
    """
    mapping = {}
    n = len(barcodes)
    # 验证specs与samples长度一致
    if len(specs) != len(samples):
        print(f"ERROR: 提供的范围数量({len(specs)})与样本名称数量({len(samples)})不一致", file=sys.stderr)
        sys.exit(1)

    # 遍历每个索引规范，将对应barcode映射到样本
    for spec, sample in zip(specs, samples):
        indices = parse_index_spec(spec, n)
        for idx in indices:
            if idx < 0 or idx >= n:
                print(f"ERROR: 索引 {idx + 1} 超出barcode列表长度{n}", file=sys.stderr)
                sys.exit(1)
            mapping[barcodes[idx]] = sample
    return mapping


def assign_sample_id(df, mapping):
    """
    针对DataFrame的CB列，提取第二段barcode并查表映射SampleID
    """
    cb = df['CB'].astype(str)
    second = cb.str.split('_', expand=True)[1]
    df['SampleID'] = second.map(mapping).fillna('Unknown')
    return df


def main():
    args = parse_args()

    # 读取barcode列表
    barcodes = load_barcodes(args.barcodelist)

    # 构建映射表
    mapping = build_mapping(barcodes, args.splitCB, args.splitSample)

    # 读取输入表格
    try:
        df = pd.read_csv(args.inputfile, sep='\t', dtype=str)
    except Exception as e:
        print(f"ERROR: 读取文件失败 {args.inputfile}: {e}", file=sys.stderr)
        sys.exit(1)

    if 'CB' not in df.columns:
        print("ERROR: 未检测到CB列", file=sys.stderr)
        sys.exit(1)

    # 添加SampleID列
    df = assign_sample_id(df, mapping)

    # 写出结果
    try:
        df.to_csv(args.output, sep='\t', index=False)
        print(f"Success: 输出文件 {args.output}")
    except Exception as e:
        print(f"ERROR: 写入文件失败 {args.output}: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()