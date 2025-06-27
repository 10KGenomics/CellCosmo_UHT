import os
import argparse
import gzip
from glob import glob
import shutil

# 解析命令行参数
def parse_args():
    parser = argparse.ArgumentParser(description="Modify barcodes.tsv.gz files based on the given flag.")
    # 是否修改 barcodes.tsv.gz 文件内容的标志
    parser.add_argument('--ModifyBarcodeName', type=lambda x: (str(x).lower() == 'true'), default=False,
                        help="Whether to modify the barcodes.tsv.gz files. Set to True or False.")
    # 输入矩阵所在的目录
    parser.add_argument('--inputMatrixDir', required=True,
                        help="Directory containing *_feature_bc_matrix folders.")
    # 输出矩阵保存的目录
    parser.add_argument('--outputMatrixDir', required=True,
                        help="Directory to save the modified matrix folders.")
    return parser.parse_args()

# 修改 barcodes.tsv.gz 文件内容的函数
def modify_barcode_file(input_file, output_file, prefix):
    with gzip.open(input_file, 'rt') as infile, gzip.open(output_file, 'wt') as outfile:
        for line in infile:
            # 在每行开头添加前缀
            new_line = f"{prefix}{line}"
            outfile.write(new_line)

# 主函数，处理整个流程
def main():
    args = parse_args()
    # 获取输入目录下所有的 *_feature_bc_matrix 文件夹
    input_lib_dirs = glob(os.path.join(args.inputMatrixDir, "*_feature_bc_matrix"))
    # 如果输入目录下没有找到符合条件的文件夹，抛出错误
    if not input_lib_dirs:
        raise FileNotFoundError(f"No *_feature_bc_matrix folders found in {args.inputMatrixDir}.")
    # 如果输出目录不存在，创建该目录
    if not os.path.exists(args.outputMatrixDir):
        os.makedirs(args.outputMatrixDir)

    for input_lib_dir in input_lib_dirs:
        # 获取当前文件夹的名称
        lib_name = os.path.basename(input_lib_dir)
        # 构建输出文件夹的路径
        output_lib_dir = os.path.join(args.outputMatrixDir, lib_name)
        # 如果输出文件夹不存在，创建该文件夹
        if not os.path.exists(output_lib_dir):
            os.makedirs(output_lib_dir)
        # 获取文件夹前缀，以 '_' 为分隔符取第一部分
        prefix = lib_name.split('_')[0] + '-'

        # 处理 barcodes.tsv.gz 文件
        barcode_input_file = os.path.join(input_lib_dir, "barcodes.tsv.gz")
        barcode_output_file = os.path.join(output_lib_dir, "barcodes.tsv.gz")
        if args.ModifyBarcodeName:
            # 如果需要修改，调用修改函数
            modify_barcode_file(barcode_input_file, barcode_output_file, prefix)
        else:
            # 如果不需要修改，直接复制文件
            shutil.copy2(barcode_input_file, barcode_output_file)

        # 复制 features.tsv.gz 和 matrix.mtx.gz 文件，这两个文件不做改动
        for file_name in ["features.tsv.gz", "matrix.mtx.gz"]:
            input_file = os.path.join(input_lib_dir, file_name)
            output_file = os.path.join(output_lib_dir, file_name)
            shutil.copy2(input_file, output_file)

        print(f"Processed {lib_name} and saved to {output_lib_dir}.")

    print("All processing completed.")

if __name__ == "__main__":
    main()
    