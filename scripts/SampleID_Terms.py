import pandas as pd
import sys

def main():
    """
    主处理函数，执行以下操作：
    1. 解析命令行参数
    2. 读取并校验输入文件
    3. 按样本分组计算统计量
    4. 计算各指标的样本间百分比
    5. 生成并保存结果文件
    """
    # 解析命令行参数
    if len(sys.argv) < 3:
        print("用法: python script.py <分组列> <统计列1> <统计列2> ...")
        sys.exit(1)
    
    group_column = sys.argv[1]    # 分组列（如SampleID）
    terms = sys.argv[2:]          # 需要统计的指标列

    # 读取输入文件
    input_file = "2-Cell_Summary_withSampleID.xls"
    try:
        df = pd.read_csv(input_file, sep="\t")
    except FileNotFoundError:
        print(f"错误：输入文件 {input_file} 未找到")
        sys.exit(1)
    except Exception as e:
        print(f"读取文件时发生错误：{str(e)}")
        sys.exit(1)

    # 校验列名有效性
    missing_cols = [col for col in [group_column] + terms if col not in df.columns]
    if missing_cols:
        print(f"错误：以下列不存在于输入文件中：{', '.join(missing_cols)}")
        sys.exit(1)

    # 按样本分组计算总和
    sum_df = df.groupby(group_column)[terms].sum().reset_index()

    # 计算各指标的样本间百分比
    percent_df = sum_df.copy()
    for term in terms:
        total = percent_df[term].sum()
        percent_df[term] = (percent_df[term] / total * 100).round(2)
    percent_df = percent_df.rename(columns={term: f"{term}_%" for term in terms})

    # 合并结果
    merged_df = pd.merge(sum_df, percent_df, on=group_column)

    # 排列列顺序（原始列后紧跟对应的百分比列）
    column_order = [group_column]
    for term in terms:
        column_order.extend([term, f"{term}_%"])
    
    result_df = merged_df[column_order]

    # 保存结果文件
    output_file = "SampleID_Summary.xls"
    result_df.to_csv(output_file, sep="\t", index=False)
    print(f"处理完成！结果已保存至：{output_file}")

if __name__ == "__main__":
    main()