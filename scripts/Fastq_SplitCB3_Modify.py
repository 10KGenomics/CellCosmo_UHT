# -*- coding: utf-8 -*-

import gzip
import argparse
import os
from collections import defaultdict
import sys
import time
import datetime
from itertools import zip_longest  # 导入zip_longest

def parse_range(range_str):
    """解析范围字符串如'1-18'，返回起始和结束索引（从1开始计数）"""
    try:
        parts = range_str.split('-')
        if len(parts) != 2:
            raise ValueError("范围格式应为'start-end'")
            
        start, end = map(int, parts)
        if start <= 0 or end <= 0:
            raise ValueError("范围值必须为正整数")
        if start > end:
            raise ValueError("起始位置不能大于结束位置")
            
        return start, end
    except Exception as e:
        raise ValueError(f"无效的范围格式'{range_str}': {e}")

def read_fastq_batch(file, batch_size=10000):
    """批量读取FASTQ记录，确保正确处理文件结束"""
    while True:
        batch = []
        records = []
        
        # 读取4*batch_size行
        for _ in range(4 * batch_size):
            line = file.readline()
            if not line:  # 文件结束
                break
            records.append(line.strip())
        
        if not records:  # 没有更多数据
            if batch:
                yield batch
            return
        
        # 分组为4行一组
        for i in range(0, len(records), 4):
            if i + 3 < len(records):
                header, seq, plus, qual = records[i:i+4]
                batch.append((header, seq, plus, qual))
            else:
                # 处理不完整的记录（警告）
                print(f"警告：发现不完整的FASTQ记录（{len(records[i:])}行）")
        
        if batch:
            yield batch

def modify_record(r1_record, r2_record, link_seq, link_start, link_end, cut_r1, cut_r2):
    """
    高效修改FASTQ记录：
    1. 减少不必要的字符串操作
    2. 优化内存分配
    3. 添加边界检查避免错误
    
    参数:
        r1_record (tuple): (header, sequence, plus, quality)
        r2_record (tuple): (header, sequence, plus, quality)
        link_seq (str): 需要匹配的link序列
        link_start (int): link序列起始位置(0-based)
        link_end (int): link序列结束位置(1-based)
        cut_r1 (int): R1保留长度
        cut_r2 (int): R2保留长度
    
    返回:
        tuple: (修改后的R1记录字符串, 修改后的R2记录字符串)
    """
    header_r1, seq_r1, plus_r1, qual_r1 = r1_record
    header_r2, seq_r2, plus_r2, qual_r2 = r2_record
    
    # 获取序列长度
    r1_len = len(seq_r1)
    r2_len = len(seq_r2)
    
    # 功能1：检查并修改R1序列
    # 仅在长度足够且序列匹配时才修改
    if r1_len >= link_end and seq_r1[link_start:link_end] == link_seq:
        # 在link_start处插入'C'，移除最后一个碱基
        # 使用一次性字符串构建避免多次拼接
        seq_modified = f"{seq_r1[:link_start]}C{seq_r1[link_start:-1]}"
        # 质量值处理：复制插入点的质量值
        qual_modified = f"{qual_r1[:link_start]}{qual_r1[link_start]}{qual_r1[link_start:-1]}"
    else:
        # 不匹配时不修改
        seq_modified, qual_modified = seq_r1, qual_r1
    
    # 功能2：序列剪切与重组
    # 边界检查避免错误
    cut_r1 = min(cut_r1, len(seq_modified))
    cut_r2 = min(cut_r2, r2_len)
    
    # 新R1 = R1前cut_r1 + R2前cut_r2
    new_r1_seq = f"{seq_modified[:cut_r1]}{seq_r2[:cut_r2]}"
    new_r1_qual = f"{qual_modified[:cut_r1]}{qual_r2[:cut_r2]}"
    
    # 新R2 = R2剩余部分
    new_r2_seq = seq_r2[cut_r2:]
    new_r2_qual = qual_r2[cut_r2:]
    
    # 返回修改后的记录字符串
    return (f"{header_r1}\n{new_r1_seq}\n{plus_r1}\n{new_r1_qual}\n",
            f"{header_r2}\n{new_r2_seq}\n{plus_r2}\n{new_r2_qual}\n")

def find_best_match(cb_barcode, barcode_set, max_mismatch):
    """优化条形码匹配算法：
    1. 使用提前终止策略减少不必要的比较
    2. 避免计算完整汉明距离
    3. 使用局部变量加速访问
    
    参数:
        cb_barcode (str): 从read中提取的条形码
        barcode_set (set): 参考条形码集合
        max_mismatch (int): 最大允许错配数
    
    返回:
        tuple: (最佳匹配的条形码, 错配数)
    """
    best_match = None
    min_mismatch = float('inf')
    barcode_len = len(cb_barcode)
    
    # 使用局部变量加速循环
    for ref_barcode in barcode_set:
        ref_len = len(ref_barcode)
        # 长度不匹配直接跳过
        if ref_len != barcode_len:
            continue
            
        # 计算错配数时提前终止
        mismatch = 0
        for i in range(barcode_len):
            if cb_barcode[i] != ref_barcode[i]:
                mismatch += 1
                # 提前终止：当错配数超过当前最小值
                if mismatch > min_mismatch:
                    break
        else:  # 循环正常结束（没有break）
            if mismatch < min_mismatch:
                min_mismatch = mismatch
                best_match = ref_barcode
                # 找到完美匹配提前退出
                if min_mismatch == 0:
                    break
            continue
    
    # 检查是否满足最大错配要求
    if min_mismatch > max_mismatch:
        return None, min_mismatch
    return best_match, min_mismatch

def format_timestamp():
    """返回格式化的当前时间戳"""
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def main():
    # 记录脚本开始时间
    start_time = time.time()
    start_timestamp = format_timestamp()
    
    # 创建参数解析器
    parser = argparse.ArgumentParser(
        description='FASTQ处理器：根据CB3条形码拆分FASTQ文件并进行序列修改',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
示例用法:
  python Fastq_SplitCB3_Modify.py \\
    --Sample WYX1 \\
    --CB3list ./CB3.ref.list \\
    --CB3_Num 1-18 \\
    --Mismatch 2 \\
    --Location 1-8 \\
    --link_location 10-15 \\
    --link_seq CAGAGC \\
    --cut_R1_length 44 \\
    --cut_R2_length 26 \\
    --inputR1 01-Fastq_Modify_ParaFly/WYX1_R1.fastq.gz \\
    --inputR2 01-Fastq_Modify_ParaFly/WYX1_R2.fastq.gz \\
    --outDir ./01-Fastq_Modify_ParaFly/WYX1
        
输出文件:
  在输出目录中为每个条形码创建两个文件:
    <Sample>_<Barcode>_R1.fastq.gz
    <Sample>_<Barcode>_R2.fastq.gz
        """
    )
    
    # 样本和条形码参数
    sample_group = parser.add_argument_group('样本和条形码参数')
    sample_group.add_argument('--Sample', required=True, 
                             help='样本名称，用于输出文件前缀（例如: WYX1）')
    sample_group.add_argument('--CB3list', required=True, 
                             help='CB3条形码列表文件路径，每行一个条形码序列')
    sample_group.add_argument('--CB3_Num', required=True, 
                             help='使用的条形码范围，格式为"start-end"（例如: 1-18表示使用前18个条形码）')
    sample_group.add_argument('--Mismatch', type=int, required=True, choices=range(0,5), 
                             help='最大允许的错配数（0-4）')
    sample_group.add_argument('--Location', required=True, 
                             help='条形码在R2 read中的位置，格式为"start-end"（例如: 1-8表示R2的第1到8个碱基）')
    
    # 序列修改参数
    modify_group = parser.add_argument_group('序列修改参数')
    modify_group.add_argument('--link_location', required=True, 
                             help='Link序列在R1 read中的位置，格式为"start-end"（例如: 10-15）')
    modify_group.add_argument('--link_seq', required=True, 
                             help='需要匹配的Link序列（例如: CAGAGC）')
    modify_group.add_argument('--cut_R1_length', type=int, required=True, 
                             help='从R1 read开头保留的碱基数')
    modify_group.add_argument('--cut_R2_length', type=int, required=True, 
                             help='从R2 read开头保留的碱基数')
    
    # 输入输出参数
    io_group = parser.add_argument_group('输入输出参数')
    io_group.add_argument('--inputR1', required=True, 
                         help='输入的R1 FASTQ.gz文件路径')
    io_group.add_argument('--inputR2', required=True, 
                         help='输入的R2 FASTQ.gz文件路径')
    io_group.add_argument('--outDir', required=True, 
                         help='输出目录路径，将在此目录中创建输出文件')
    
    # 性能参数
    perf_group = parser.add_argument_group('性能参数')
    perf_group.add_argument('--batch_size', type=int, default=10000, 
                           help='批量处理大小（建议10000-50000，较大的值可提高速度但增加内存使用）')
    perf_group.add_argument('--compress_level', type=int, default=1, choices=range(1,10),
                           help='输出文件的压缩级别（1=最快压缩，9=最高压缩）')
    perf_group.add_argument('--progress_interval', type=int, default=1000000, 
                           help='进度报告间隔（处理的read数）')
    args = parser.parse_args()
    
    print(f"[{format_timestamp()}] 开始处理样本: {args.Sample}")
    print(f"[{format_timestamp()}] 参数设置:")
    print(f"  样本名称: {args.Sample}")
    print(f"  CB3列表文件: {args.CB3list}")
    print(f"  CB3范围: {args.CB3_Num}")
    print(f"  最大错配: {args.Mismatch}")
    print(f"  条形码位置: {args.Location}")
    print(f"  Link序列位置: {args.link_location}")
    print(f"  Link序列: {args.link_seq}")
    print(f"  R1保留长度: {args.cut_R1_length}")
    print(f"  R2保留长度: {args.cut_R2_length}")
    print(f"  输入R1文件: {args.inputR1}")
    print(f"  输入R2文件: {args.inputR2}")
    print(f"  输出目录: {args.outDir}")
    print(f"  批处理大小: {args.batch_size}")
    print(f"  压缩级别: {args.compress_level}")
    print(f"  进度报告间隔: {args.progress_interval} reads")

    # 创建输出目录
    os.makedirs(args.outDir, exist_ok=True)
    print(f"[{format_timestamp()}] 创建输出目录: {args.outDir}")

    # 解析条形码范围
    try:
        cb_start, cb_end = parse_range(args.CB3_Num)
        print(f"[{format_timestamp()}] 解析CB3范围: {cb_start}-{cb_end}")
    except ValueError as e:
        print(f"[{format_timestamp()}] 错误: {e}")
        sys.exit(1)
    
    # 解析条形码位置
    try:
        loc_start, loc_end = parse_range(args.Location)
        start_idx = loc_start - 1  # 转换为0-based索引
        end_idx = loc_end
        print(f"[{format_timestamp()}] 解析条形码位置: {loc_start}-{loc_end} (0-based索引: {start_idx}-{end_idx})")
    except ValueError as e:
        print(f"[{format_timestamp()}] 错误: {e}")
        sys.exit(1)
    
    # 解析Link序列位置
    try:
        link_start, link_end = parse_range(args.link_location)
        link_start0 = link_start - 1  # 转换为0-based索引
        print(f"[{format_timestamp()}] 解析Link位置: {link_start}-{link_end} (0-based索引: {link_start0}-{link_end})")
        
        # 验证Link序列长度
        if (link_end - link_start) != len(args.link_seq):
            print(f"[{format_timestamp()}] 警告: Link序列长度({len(args.link_seq)})与位置范围({link_end-link_start})不匹配")
    except ValueError as e:
        print(f"[{format_timestamp()}] 错误: {e}")
        sys.exit(1)

    # 读取CB3条形码
    selected_barcodes = []
    barcode_lengths = set()
    print(f"[{format_timestamp()}] 从 '{args.CB3list}' 读取条形码...")
    
    try:
        with open(args.CB3list, 'r') as f:
            for i, line in enumerate(f, 1):
                if cb_start <= i <= cb_end:
                    barcode = line.strip().upper()
                    if not barcode:  # 跳过空行
                        continue
                    selected_barcodes.append(barcode)
                    barcode_lengths.add(len(barcode))
    except Exception as e:
        print(f"[{format_timestamp()}] 读取条形码文件错误: {e}")
        sys.exit(1)
    
    if not selected_barcodes:
        print(f"[{format_timestamp()}] 错误: 在指定范围 {cb_start}-{cb_end} 内未找到条形码")
        sys.exit(1)
    
    # 检查条形码长度一致性
    if len(barcode_lengths) > 1:
        print(f"[{format_timestamp()}] 警告: 条形码长度不一致 {barcode_lengths}，可能导致匹配问题")
    
    print(f"[{format_timestamp()}] 从库中选择了 {len(selected_barcodes)} 个条形码")

    # 为每个条形码创建输出文件
    output_files = {}
    print(f"[{format_timestamp()}] 为每个条形码创建输出文件...")
    
    try:
        for barcode in selected_barcodes:
            base_name = os.path.join(args.outDir, f"{args.Sample}_{barcode}")
            r1_path = f"{base_name}_R1.fastq.gz"
            r2_path = f"{base_name}_R2.fastq.gz"
            
            r1_out = gzip.open(r1_path, 'wt', compresslevel=args.compress_level)
            r2_out = gzip.open(r2_path, 'wt', compresslevel=args.compress_level)
            output_files[barcode] = (r1_out, r2_out)
            
            print(f"  创建: {r1_path} 和 {r2_path}")
    except Exception as e:
        print(f"[{format_timestamp()}] 创建输出文件错误: {e}")
        sys.exit(1)

    # 初始化统计信息
    stats = defaultdict(int)
    total_reads = 0
    matched_reads = 0
    last_report = 0
    batch_count = 0

    # 关键优化：将条形码列表转换为集合（更快查找）
    barcode_set = set(selected_barcodes)
    
    print(f"[{format_timestamp()}] 开始处理FASTQ文件...")
    
    try:
        with gzip.open(args.inputR1, 'rt') as f_r1, gzip.open(args.inputR2, 'rt') as f_r2:
            # 创建批量读取器
            r1_batches = read_fastq_batch(f_r1, args.batch_size)
            r2_batches = read_fastq_batch(f_r2, args.batch_size)
            
            # 使用zip_longest处理不同长度的输入
            for batch_idx, batch_pair in enumerate(zip_longest(r1_batches, r2_batches, fillvalue=None)):
                r1_batch, r2_batch = batch_pair
                
                # 检查文件是否结束
                if r1_batch is None or r2_batch is None:
                    if r1_batch is None and r2_batch is not None:
                        print(f"[{format_timestamp()}] 错误: R1文件结束但R2还有数据 (批次 {batch_idx})")
                    elif r2_batch is None and r1_batch is not None:
                        print(f"[{format_timestamp()}] 错误: R2文件结束但R1还有数据 (批次 {batch_idx})")
                    else:
                        # 两个文件同时结束
                        print(f"[{format_timestamp()}] 文件处理完成 (批次 {batch_idx})")
                    break
                
                batch_count += 1
                batch_size = len(r1_batch)
                total_reads += batch_size
                
                # 进度报告
                if total_reads - last_report >= args.progress_interval:
                    elapsed = time.time() - start_time
                    reads_per_sec = total_reads / elapsed if elapsed > 0 else 0
                    print(f"[{format_timestamp()}] 已处理: {total_reads:,} reads | "
                          f"速度: {reads_per_sec:,.0f} reads/sec | "
                          f"批次: {batch_count}")
                    last_report = total_reads
                    
                    # 刷新输出缓冲区
                    for r1_out, r2_out in output_files.values():
                        r1_out.flush()
                        r2_out.flush()
                
                # 为每个条形码初始化本批次缓冲区
                batch_buffers = {barcode: ([], []) for barcode in barcode_set}
                
                # 处理本批次中的每条记录
                for r1_record, r2_record in zip(r1_batch, r2_batch):
                    # 从R2提取条形码
                    r2_seq = r2_record[1]
                    r2_len = len(r2_seq)
                    
                    # 检查R2长度是否足够
                    if r2_len < end_idx:
                        stats['too_short'] += 1
                        continue
                    
                    # 提取条形码区域
                    cb_barcode = r2_seq[start_idx:end_idx].upper()
                    
                    # 使用优化的匹配函数
                    best_match, min_mismatch = find_best_match(
                        cb_barcode, 
                        barcode_set, 
                        args.Mismatch
                    )
                    
                    # 检查是否匹配
                    if best_match is not None:
                        matched_reads += 1
                        stats[best_match] += 1
                        
                        # 修改序列
                        mod_r1, mod_r2 = modify_record(
                            r1_record, r2_record,
                            args.link_seq, link_start0, link_end,
                            args.cut_R1_length, args.cut_R2_length
                        )
                        
                        # 添加到缓冲区
                        r1_buffer, r2_buffer = batch_buffers[best_match]
                        r1_buffer.append(mod_r1)
                        r2_buffer.append(mod_r2)
                    else:
                        stats['unmatched'] += 1
                
                # 将本批次数据写入文件
                for barcode in barcode_set:
                    r1_list, r2_list = batch_buffers[barcode]
                    if r1_list:
                        r1_out, r2_out = output_files[barcode]
                        # 使用join一次性写入提高效率
                        r1_out.write(''.join(r1_list))
                        r2_out.write(''.join(r2_list))
    except Exception as e:
        print(f"[{format_timestamp()}] 处理过程中发生错误: {e}")
        # 确保在异常时关闭所有文件
        for barcode, (r1_out, r2_out) in output_files.items():
            try:
                r1_out.close()
                r2_out.close()
            except:
                pass
        sys.exit(1)
    
    finally:  # 确保在任何情况下都关闭文件
        print(f"[{format_timestamp()}] 确保关闭所有输出文件...")
        for barcode, (r1_out, r2_out) in output_files.items():
            try:
                r1_out.close()
                print(f"  已关闭: {args.Sample}_{barcode}_R1.fastq.gz")
            except Exception as e:
                print(f"  关闭 {barcode} R1 文件时出错: {e}")
            try:
                r2_out.close()
                print(f"  已关闭: {args.Sample}_{barcode}_R2.fastq.gz")
            except Exception as e:
                print(f"  关闭 {barcode} R2 文件时出错: {e}")
    
    # 计算总运行时间
    end_time = time.time()
    total_time = end_time - start_time
    hours, remainder = divmod(total_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    
    # 打印统计信息
    print(f"\n[{format_timestamp()}] 处理完成!")
    print("=" * 80)
    print(f"总运行时间: {int(hours)}小时 {int(minutes)}分钟 {seconds:.2f}秒")
    print(f"总处理reads数: {total_reads:,}")
    
    if total_reads > 0:
        reads_per_sec = total_reads / total_time
        match_percent = matched_reads / total_reads * 100
        print(f"处理速度: {reads_per_sec:,.0f} reads/sec")
        print(f"匹配到条形码的reads: {matched_reads:,} ({match_percent:.2f}%)")
    else:
        print("警告: 未处理任何reads")
    
    print("\n条形码分布:")
    for barcode in selected_barcodes:
        count = stats.get(barcode, 0)
        percent = count / matched_reads * 100 if matched_reads > 0 else 0
        print(f"  {barcode}: {count:,} reads ({percent:.2f}%)")
    
    print("\n其他统计:")
    print(f"  未匹配的reads: {stats.get('unmatched', 0):,}")
    print(f"  长度不足的reads: {stats.get('too_short', 0):,}")
    print("=" * 80)

    # 明确退出
    sys.exit(0)

if __name__ == "__main__":
    main()