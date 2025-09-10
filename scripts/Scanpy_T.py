import os
from glob import glob
import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import re
import logging
from collections import defaultdict
from typing import Dict, Set, Tuple, List  # 添加 List 导入

# 配置日志系统
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger(__name__)

##########################################
# 计算 AnnData 对象的基本统计指标
##########################################
def compute_library_stats(adata):
    """
    计算单个 AnnData 对象的统计指标：
      - Num_Cells: 细胞数
      - Total_UMI: 全部转录本（UMI）数
      - Median_UMI_per_Cell: 每个细胞的 UMI 中位数
      - Mean_UMI_per_Cell: 每个细胞的 UMI 平均数
      - Median_Genes_per_Cell: 每个细胞检测到的基因数中位值
      - Mean_Genes_per_Cell: 每个细胞检测到的基因数均值
      - Total_Genes_Covered: 文库中总共检测到的基因种类数
    """
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    num_cells = adata.n_obs
    total_umi = adata.obs['total_counts'].sum()
    median_umi = np.median(adata.obs['total_counts'])
    mean_umi = np.mean(adata.obs['total_counts'])
    median_genes = np.median(adata.obs['n_genes_by_counts'])
    mean_genes = np.mean(adata.obs['n_genes_by_counts'])

    if hasattr(adata.X, 'toarray'):
        gene_totals = np.array(adata.X.toarray().sum(axis=0)).flatten()
    else:
        gene_totals = np.array(adata.X.sum(axis=0)).flatten()
    gene_covered = np.sum(gene_totals > 0)
    return {
        'Num_Cells': num_cells,
        'Total_UMI': total_umi,
        'Median_UMI_per_Cell': median_umi,
        'Mean_UMI_per_Cell': mean_umi,
        'Median_Genes_per_Cell': median_genes,
        'Mean_Genes_per_Cell': mean_genes,
        'Total_Genes_Covered': gene_covered
    }

def parse_args():
    parser = argparse.ArgumentParser(
        description="单细胞矩阵处理流程：支持raw/filtered矩阵的读取、合并和样本拆分"
    )
    parser.add_argument(
        '--inputMatrixDir', '-i', required=True,
        help="10X矩阵父目录，包含raw和filtered子目录，例如：04-Matrix/"
    )
    parser.add_argument(
        '--matrixType', '-m', required=True, choices=['raw', 'filtered'],
        help="矩阵类型选择：raw 或 filtered"
    )
    parser.add_argument(
        '--inputfileBarcodelist', '-b', required=True,
        help="mapping.txt 文件路径（包含barcode和文库信息）"
    )
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
        help="与splitCB对应的样本名称列表（顺序必须一致）"
    )
    parser.add_argument(
        '--outputDir', '-o', default='.',
        help="输出目录路径，默认为当前目录"
    )
    return parser.parse_args()

def parse_index_spec(spec: str) -> Set[int]:
    """
    解析索引范围规范，返回1-based索引集合
    
    支持格式:
    - 单个数字: "1"
    - 范围: "1-64"
    - 混合: "1,3,5-10"
    - 等差数列: "[1,100,5]"
    
    注意: 返回的是1-based索引值集合
    """
    indices = set()
    
    # 处理等差数列语法 [start,end,step]
    if spec.startswith('[') and spec.endswith(']'):
        try:
            parts = spec[1:-1].split(',')
            if len(parts) != 3:
                raise ValueError("需要3个参数: start,end,step")
                
            start = int(parts[0].strip())
            end = int(parts[1].strip())
            step = int(parts[2].strip())
            
            if step == 0:
                raise ValueError("步长不能为0")
                
            # 生成数列
            current = start
            if start <= end:
                while current <= end:
                    indices.add(current)
                    current += abs(step)
            else:
                while current >= end:
                    indices.add(current)
                    current -= abs(step)
                    
            return indices
            
        except Exception as e:
            raise ValueError(f"无效的等差数列语法 '{spec}': {str(e)}")
    
    # 处理常规语法 (逗号分隔的单个数字或范围)
    parts = [p.strip() for p in spec.split(',') if p.strip()]
    
    for part in parts:
        if '-' in part:
            try:
                a, b = map(int, part.split('-'))
                low, high = sorted([a, b])
                indices.update(range(low, high + 1))
            except:
                raise ValueError(f"无效范围格式 '{part}'")
        else:
            try:
                indices.add(int(part))
            except:
                raise ValueError(f"无效索引值 '{part}'")
    
    return indices

def parse_bc2_info(info_str: str) -> int:
    """
    从文库信息字符串中提取CB2索引号（1-based）
    
    示例输入: "bc20,a5,l1,bc2-135"
    返回: 135
    """
    # 查找 bc2-xxx 模式
    match = re.search(r'bc2-(\d+)', info_str)
    if match:
        return int(match.group(1))
    
    # 查找其他可能格式
    match = re.search(r'bc2_(\d+)', info_str)
    if match:
        return int(match.group(1))
    
    # 尝试更宽松的匹配
    match = re.search(r'bc2.*?(\d+)', info_str)
    if match:
        return int(match.group(1))
    
    raise ValueError(f"未找到有效的CB2索引号: {info_str}")

def create_sample_mapping(
    mapping_file: str, 
    split_specs: List[str], 
    samples: List[str]
) -> Dict[str, str]:
    """
    创建16bp barcode到样本的映射
    
    步骤:
    1. 解析每个样本的索引范围规范
    2. 读取mapping.txt文件
    3. 为每个barcode解析CB2索引号
    4. 根据索引范围分配样本
    
    返回:
    - barcode_to_sample: 16bp barcode到样本名称的映射字典
    """
    # 校验参数
    if len(split_specs) != len(samples):
        raise ValueError("splitCB和splitSample数量不匹配")
    
    # 解析每个样本的索引范围
    sample_indices = {}
    for spec, sample in zip(split_specs, samples):
        try:
            sample_indices[sample] = parse_index_spec(spec)
        except ValueError as e:
            raise ValueError(f"样本 '{sample}' 的索引规范错误: {str(e)}")
    
    # 创建barcode到样本的映射
    barcode_to_sample = {}
    bc2_index_count = defaultdict(int)
    unmapped_barcodes = 0
    
    try:
        with open(mapping_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    logger.warning(f"跳过无效行 {line_num}: 列数不足")
                    continue
                
                barcode = parts[0].strip()
                info_str = parts[1].strip()
                
                try:
                    # 解析CB2索引号
                    bc2_index = parse_bc2_info(info_str)
                    bc2_index_count[bc2_index] += 1
                    
                    # 查找对应的样本
                    found = False
                    for sample, indices in sample_indices.items():
                        if bc2_index in indices:
                            barcode_to_sample[barcode] = sample
                            found = True
                            break
                    
                    if not found:
                        logger.debug(f"CB2索引 {bc2_index} 不属于任何样本范围")
                        unmapped_barcodes += 1
                
                except ValueError as e:
                    logger.warning(f"行 {line_num} 解析错误: {str(e)}")
                    unmapped_barcodes += 1
    
    except IOError as e:
        raise IOError(f"无法读取mapping文件: {str(e)}")
    
    # 检查是否有barcode被分配
    if not barcode_to_sample:
        raise RuntimeError("没有barcode被分配到样本，请检查mapping.txt和splitCB参数")
    
    # 记录统计信息
    total_mapped = len(barcode_to_sample)
    logger.info(f"成功映射 {total_mapped} 个barcode到样本")
    logger.info(f"未映射的barcode数量: {unmapped_barcodes}")
    
    return barcode_to_sample

def main():
    args = parse_args()

    # 参数校验
    if len(args.splitCB) != len(args.splitSample):
        raise ValueError(f"splitCB数量({len(args.splitCB)})与splitSample数量({len(args.splitSample)})不匹配")
    
    # 创建输出目录结构
    combined_dir = os.path.join(args.outputDir, f"05-Combined/{args.matrixType}")
    sample_dir = os.path.join(args.outputDir, f"06-CB2_Sample/{args.matrixType}")
    os.makedirs(combined_dir, exist_ok=True)
    os.makedirs(sample_dir, exist_ok=True)

    # 步骤1: 创建barcode到样本的映射
    logger.info("创建barcode到样本的映射...")
    try:
        barcode_to_sample = create_sample_mapping(
            args.inputfileBarcodelist,
            args.splitCB,
            args.splitSample
        )
    except Exception as e:
        logger.error(f"映射创建失败: {str(e)}")
        sys.exit(1)
    
    # 根据矩阵类型确定目录后缀
    matrix_suffix = {
        'raw': 'raw_feature_bc_matrix',
        'filtered': 'filtered_feature_bc_matrix'
    }[args.matrixType]
    
    # 获取矩阵目录列表
    matrix_dir = os.path.join(args.inputMatrixDir)
    lib_dirs = glob(os.path.join(matrix_dir, f"*_{matrix_suffix}"))
    if not lib_dirs:
        raise FileNotFoundError(f"未找到{args.matrixType}矩阵目录: {matrix_dir}/*_{matrix_suffix}")
    
    logger.info(f"找到 {len(lib_dirs)} 个{args.matrixType}矩阵目录")

    # 处理每个文库
    metadata_list = []
    adata_list = []
    for lib_dir in lib_dirs:
        lib_name = os.path.basename(lib_dir).replace(f"_{matrix_suffix}", "")
        logger.info(f"处理文库: {lib_name} ({args.matrixType}矩阵)")
        
        # 读取10x矩阵
        try:
            adata = sc.read_10x_mtx(
                lib_dir, 
                var_names='gene_symbols', 
                cache=True,
                make_unique=True  # 确保基因名唯一
            )
        except Exception as e:
            logger.error(f"读取文库 {lib_name} 失败: {str(e)}")
            continue
        
        # 添加文库信息
        adata.obs['Library'] = lib_name
        
        # 提取16bp barcode（去掉前缀）
        try:
            # 分割barcode名称，获取第二部分
            adata.obs['16bp_Barcode'] = adata.obs_names.str.split('-').str[1]
        except Exception as e:
            logger.error(f"处理文库 {lib_name} 的barcode失败: {str(e)}")
            logger.error(f"示例barcode: {adata.obs_names[0] if adata.n_obs > 0 else 'N/A'}")
            continue
        
        # 分配样本
        adata.obs['Sample'] = adata.obs['16bp_Barcode'].map(barcode_to_sample)
        
        # 检查未匹配的barcode
        unmatched = adata.obs['Sample'].isnull().sum()
        if unmatched > 0:
            logger.warning(f"文库 {lib_name} 中有 {unmatched} 个barcode未匹配到样本")
        
        # 计算统计信息
        try:
            lib_stats = compute_library_stats(adata)
            lib_stats['Library'] = lib_name
            lib_stats['Unmatched_Barcodes'] = unmatched
            metadata_list.append(lib_stats)
        except Exception as e:
            logger.error(f"计算文库 {lib_name} 统计信息失败: {str(e)}")
            continue
        
        adata_list.append(adata)
    
    # 如果没有成功读取任何文库，退出
    if len(adata_list) == 0:
        logger.error("没有成功读取任何文库数据，程序终止")
        sys.exit(1)
    
    # 保存单文库元数据
    perlib_meta_path = os.path.join(combined_dir, f"3-PerLib_Metadata_{args.matrixType}.xlsx")
    pd.DataFrame(metadata_list).to_excel(perlib_meta_path, index=False)
    logger.info(f"单文库元数据已保存: {perlib_meta_path}")

    # 合并所有文库
    logger.info("合并所有文库...")
    try:
        merged_adata = sc.concat(
            adata_list, 
            join='outer', 
            label='Library',
            keys=[d.obs['Library'].iloc[0] for d in adata_list],
            index_unique=None # '-' barcode添加后缀
        )
    except Exception as e:
        logger.error(f"合并文库失败: {str(e)}")
        sys.exit(1)
    
    # 计算并保存合并后元数据
    merge_meta_path = os.path.join(combined_dir, f"3-Merge_Metadata_{args.matrixType}.xlsx")
    try:
        merge_stats = compute_library_stats(merged_adata)
        merge_stats['Total_Libraries'] = len(lib_dirs)
        pd.DataFrame([merge_stats]).to_excel(merge_meta_path, index=False)
    except Exception as e:
        logger.error(f"计算合并后元数据失败: {str(e)}")
    
    # 保存合并的h5ad文件
    merged_h5ad = os.path.join(combined_dir, f"Merge_{args.matrixType}.h5ad")
    try:
        merged_adata.write(merged_h5ad)
        logger.info(f"合并文库已保存: {merged_h5ad}")
    except Exception as e:
        logger.error(f"保存合并文库失败: {str(e)}")

    # 按样本拆分并保存
    logger.info("按样本拆分数据...")
    sample_metadata_list = []
    
    for sample in args.splitSample:
        # 提取样本数据
        sample_mask = merged_adata.obs['Sample'] == sample
        if sample_mask.sum() == 0:
            logger.warning(f"样本 {sample} 中没有细胞，跳过")
            continue
            
        sub_adata = merged_adata[sample_mask].copy()
        
        # 保存样本数据
        sample_h5ad = os.path.join(sample_dir, f"{sample}_{args.matrixType}.h5ad")
        try:
            sub_adata.write(sample_h5ad)
        except Exception as e:
            logger.error(f"保存样本 {sample} 数据失败: {str(e)}")
            continue
        
        # 计算并保存样本元数据
        try:
            sub_metadata = compute_library_stats(sub_adata)
            sub_metadata['Sample'] = sample
            sample_meta_path = os.path.join(sample_dir, f"3-{sample}_Metadata_{args.matrixType}.xlsx")
            pd.DataFrame([sub_metadata]).to_excel(sample_meta_path, index=False)
            
            sample_metadata_list.append(sub_metadata)
            logger.info(f"样本 {sample} 数据已保存: {sample_h5ad}")
        except Exception as e:
            logger.error(f"计算样本 {sample} 统计信息失败: {str(e)}")

    # 保存所有样本的汇总元数据
    if sample_metadata_list:
        sample_summary_path = os.path.join(sample_dir, f"3-Sample_Summary_{args.matrixType}.xlsx")
        pd.DataFrame(sample_metadata_list).to_excel(sample_summary_path, index=False)
    
    # 最终报告
    logger.info("=" * 60)
    logger.info(f"处理完成! 矩阵类型: {args.matrixType}")
    logger.info(f"总处理文库数: {len(lib_dirs)}")
    logger.info(f"总细胞数: {merged_adata.n_obs}")
    logger.info(f"输出目录: {combined_dir} 和 {sample_dir}")
    logger.info("=" * 60)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logger.exception("处理过程中发生致命错误")
        sys.exit(1)