import argparse
import gzip
from itertools import islice
from functools import partial

def read_fastq_batch(file, batch_size=1000):
    """批量读取 FASTQ 记录，减少 I/O 调用次数"""
    while True:
        batch = []
        # 一次性读取 4 * batch_size 行
        lines = list(islice(file, 4 * batch_size))
        if not lines:
            break
        # 按每4行分割为记录
        for i in range(0, len(lines), 4):
            header = lines[i].strip()
            seq = lines[i+1].strip()
            plus = lines[i+2].strip()
            qual = lines[i+3].strip()
            batch.append((header, seq, plus, qual))
        yield batch

def process_batch(r1_batch, r2_batch, link_seq, link_start, link_end, cut_r1, cut_r2):
    """批量处理记录，使用动态的 link_start/link_end 提升灵活性"""
    processed_r1, processed_r2 = [], []
    for r1, r2 in zip(r1_batch, r2_batch):
        header_r1, seq_r1, plus_r1, qual_r1 = r1
        header_r2, seq_r2, plus_r2, qual_r2 = r2

        # 优化点1：根据用户提供的 link_start/link_end 动态检查序列
        if len(seq_r1) >= link_end and seq_r1[link_start:link_end] == link_seq:
            # 优化点2：插入位置与 link_start 关联，而非硬编码
            # 确保序列长度足够进行操作
            if len(seq_r1) >= link_start + 1:
                # 在 link_start 位置插入 'C'，并移除最后一个碱基
                seq_modified = seq_r1[:link_start] + 'C' + seq_r1[link_start:-1]
                # 质量值处理：取 link_start 位置的质量值并保留后续（移除最后一位）
                qual_modified = qual_r1[:link_start] + qual_r1[link_start:link_start+1] + qual_r1[link_start:-1]
            else:
                # 长度不足时直接移除最后一个碱基
                seq_modified = seq_r1[:-1]
                qual_modified = qual_r1[:-1]
        else:
            # 未找到 link_seq 时不修改
            seq_modified, qual_modified = seq_r1, qual_r1

        # 优化点3：使用更高效的切片操作，避免中间变量
        new_r1_seq = (seq_modified[:cut_r1] + seq_r2[:cut_r2])[:cut_r1 + cut_r2]
        new_r1_qual = (qual_modified[:cut_r1] + qual_r2[:cut_r2])[:cut_r1 + cut_r2]
        new_r2_seq = seq_r2[cut_r2:]
        new_r2_qual = qual_r2[cut_r2:]

        # 构建输出记录
        processed_r1.append(f"{header_r1}\n{new_r1_seq}\n{plus_r1}\n{new_r1_qual}\n")
        processed_r2.append(f"{header_r2}\n{new_r2_seq}\n{plus_r2}\n{new_r2_qual}\n")
    return processed_r1, processed_r2

def main():
    parser = argparse.ArgumentParser(description='Optimized FASTQ processor')
    # 必需参数
    parser.add_argument('--input_R1', required=True, help='Input R1 FASTQ file (gzipped)')
    parser.add_argument('--input_R2', required=True, help='Input R2 FASTQ file (gzipped)')
    parser.add_argument('--link_location', required=True, help='Link序列位置，格式：start-end (1-based，如10-15)')
    parser.add_argument('--link_seq', required=True, help='Link序列，如CAGAGC')
    parser.add_argument('--cut_R1_length', type=int, required=True, help='R1保留长度')
    parser.add_argument('--cut_R2_length', type=int, required=True, help='R2保留长度')
    parser.add_argument('--output_R1', required=True, help='Output R1 FASTQ file (gzipped)')
    parser.add_argument('--output_R2', required=True, help='Output R2 FASTQ file (gzipped)')
    
    # 优化点4：添加性能调优参数
    parser.add_argument('--batch_size', type=int, default=5000, 
                       help='批量处理大小，根据内存调整 (默认:5000)')
    parser.add_argument('--compress_level', type=int, default=1,
                       help='输出压缩级别 (1速度最快，9压缩率最高，默认:1)')
    
    args = parser.parse_args()

    # 解析 link_location 参数
    try:
        start_str, end_str = args.link_location.split('-')
        link_start = int(start_str) - 1  # 转换为 0-based
        link_end = int(end_str)
        # 验证位置有效性
        if link_start < 0 or link_end <= link_start:
            raise ValueError("link_location 必须满足 start < end 且 start ≥ 1 (1-based)")
        if link_end - link_start != len(args.link_seq):
            raise ValueError("link_seq 长度与 link_location 范围不匹配")
    except Exception as e:
        raise ValueError(f"解析 --link_location 失败: {e}")

    # 优化点5：使用可配置的压缩级别和批量大小
    with gzip.open(args.input_R1, 'rt') as f_r1, \
         gzip.open(args.input_R2, 'rt') as f_r2, \
         gzip.open(args.output_R1, 'wt', compresslevel=args.compress_level) as f_out_r1, \
         gzip.open(args.output_R2, 'wt', compresslevel=args.compress_level) as f_out_r2:

        # 初始化批量读取器
        r1_reader = read_fastq_batch(f_r1, args.batch_size)
        r2_reader = read_fastq_batch(f_r2, args.batch_size)

        # 批量处理
        for r1_batch, r2_batch in zip(r1_reader, r2_reader):
            processed_r1, processed_r2 = process_batch(
                r1_batch, r2_batch, 
                args.link_seq, link_start, link_end,
                args.cut_R1_length, args.cut_R2_length
            )
            # 批量写入
            f_out_r1.write(''.join(processed_r1))
            f_out_r2.write(''.join(processed_r2))

if __name__ == "__main__":
    main()
