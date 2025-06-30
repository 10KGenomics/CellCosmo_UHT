import os
import glob
import subprocess
import argparse
from argparse import RawTextHelpFormatter

def main():
    # 解析命令行参数
    parser = argparse.ArgumentParser(
        description="10K 超高通量单细胞转录组分析流程Pipeline ",
        formatter_class=RawTextHelpFormatter
    )
    # parser.add_argument("--condaName", required=True, help="conda环境名称，比如：CellCosmo_UHT")
    parser.add_argument("--script_path", required=True, help="pipeline每步脚本存放路径，比如：/path/CellCosmo_UHT")
    
    # 样本及基础参数
    parser.add_argument("--SampleName", required=True, help="下机文库名，单个或多个(多个文库名则空格隔开)都可能，比如：'WT1-1 WT1-2'")
    parser.add_argument("--STARindex", required=True, help="参考基因组STAR index 路径，比如：/path/GRCh38_index")
    parser.add_argument("--TopCells", type=int, required=True, help="CB3单个子文库的强制细胞数,根据实验信息填写，比如：2000")
    parser.add_argument("--Threads", type=int, required=True, help="STARsolo单个任务线程数，比如：16")
    parser.add_argument("--CB3_Num", required=True, help="CB3选取范围index，比如：1-20")
    
    # Fastq处理参数
    parser.add_argument("--link1_location", default='10-15', help="link1位置，默认：10-15 ，如不是默认值则需要使用该参数，其他参数也如此")
    parser.add_argument("--link1_seq", default='CAGAGC', help="link1序列，默认：CAGAGC")
    parser.add_argument("--cut_R1_length", type=int, default=44, help="R1剪切长度，默认：44")
    parser.add_argument("--cut_R2_length", type=int, default=26, help="R2剪切长度，默认：26")
    parser.add_argument("--ParallelCPU_Fastq_Modify", type=int, default=1, help="FastqModify平行任务数，默认：1")
    
    # STARsolo参数
    parser.add_argument("--ParallelCPU_STARsolo", type=int, default=1, help="STARsolo平行任务数，默认：1")
    parser.add_argument("--STARsolo_param", help="STARsolo增加参数，比如：'--outSAMtype SAM --outReadsUnmapped Fastx'")
    
    # 混样拆分参数
    parser.add_argument("--SplitLibrary", default='False', help="下机多组fastq.gz,实际1个大样本。若单个下机文库名，则填写False；若多个下机文库名，则填写批次名；默认：False")
    # parser.add_argument("--CB2_list", required=True, help="CB2 barcode列表文件绝对路径，比如：/path/CellCosmo_UHT/Barcode/CB2.list ")
    parser.add_argument("--splitCB", required=True, help="混样CB2分类规则：x-y（1-64 65-128 129-192）或[起始值,终止值,公差值]（[1,89,8] [2,90,8] [3,91,8]）或x,y,z w,a,b（1,2,3,4 5,6,7,8 9,10,11,12） （以空格分隔不同样本，逗号分隔样本内部CB2的元素）")
    parser.add_argument("--splitSample", required=True, help="混样拆分后的样本名列表,比如：'A B C'")
    parser.add_argument("--splitSample_EstimatedCell_list", required=True, help="各样本预估细胞数列表文件,与splitSample对应数目，比如：'4-EstimatedCell_A.list 4-EstimatedCell_B.list  4-EstimatedCell_C.list'")
    parser.add_argument("--splitSample_EstimatedCellMatrix", required=True, help="各样本矩阵文件前缀,与splitSample对应数目，比如：'A_filtered_feature_bc_matrix B_filtered_feature_bc_matrix C_filtered_feature_bc_matrix'")
    args = parser.parse_args()

    # 初始化环境
    subprocess.run(f"mkdir -p log/", shell=True)
    subprocess.run(f"mkdir -p ParaFly/", shell=True)
    # os.makedirs("log", exist_ok=True)
    # print("===== 环境初始化完成 =====")

    # ------------------------ 步骤1：Fastq处理 ------------------------
    print("\n===== 开始 Fastq Modify 处理 =====")
    fastq_files = glob.glob("01-Fastq_Modify_ParaFly/*.fastq.gz")
    if not fastq_files:
        cmd = [
            f"nohup bash {args.script_path}/scripts/Fastq_Modify.sh",
            args.link1_location,
            args.link1_seq,
            str(args.cut_R1_length),
            str(args.cut_R2_length),
            str(args.ParallelCPU_Fastq_Modify),
            f"{args.script_path}/scripts",
            f"> log/out-01.Fastq_Modify.log"
        ]
        subprocess.run(" ".join(cmd), shell=True)
        print("Fastq处理已启动，日志见 log/out-01.Fastq_Modify.log")
    else:
        print("检测到已存在Fastq处理结果，跳过此步骤")

    # ------------------------ 步骤2：STARsolo分析 ------------------------
    print("\n===== 开始 STARsolo 分析 =====")
    star_logs = glob.glob("02-STARsolo/*-Log.final.out")
    if not star_logs:
        if args.STARsolo_param:
            # 正确构建命令参数，不使用shell=True
            cmd = [
                "nohup", "bash", f"{args.script_path}/scripts/STARsolo.sh",
                args.STARindex,
                str(args.TopCells),
                str(args.Threads),
                args.CB3_Num,
                str(args.ParallelCPU_STARsolo),
                args.script_path,
                args.STARsolo_param
            ]
            # 处理日志重定向
            print(f"STARsolo已启动，已设置增设STARsolo参数: {args.STARsolo_param}，日志见 log/out-02.STARsolo.log")
            with open("log/out-02.STARsolo.log", "w") as log_file:
                subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT, shell=False)
        else:
            cmd = [
                "nohup", "bash", f"{args.script_path}/scripts/STARsolo.sh",
                args.STARindex,
                str(args.TopCells),
                str(args.Threads),
                args.CB3_Num,
                str(args.ParallelCPU_STARsolo),
                args.script_path
            ]
            print("STARsolo已启动，日志见 log/out-02.STARsolo.log")
            with open("log/out-02.STARsolo.log", "w") as log_file:
                subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT, shell=False)
    else:
        print("检测到已存在STARsolo结果，跳过此步骤")

    # ------------------------ 步骤3：细胞读数汇总 ------------------------
    print("\n===== 开始 细胞读数汇总 =====")
    os.makedirs("03-CellReadsSummary", exist_ok=True)
    samples = subprocess.check_output(
        "ls -1 02-STARsolo/*-Log.out | awk -F '-Log.out' '{print $1}' | awk -F '/' '{print $2}'",
        shell=True
    ).decode().split()
    
    for sample in samples:
        # 运行CellReads_Summary.py
        subprocess.run([
            "python", f"{args.script_path}/scripts/CellReads_Summary.py",
            f"02-STARsolo/{sample}-Solo.out/GeneFull_Ex50pAS/CellReads.stats",
            f"03-CellReadsSummary/{sample}-CellReads_Summary.xls",
            sample
        ])
        
        # 运行Summary_merge.py
        subprocess.run([
            "python", f"{args.script_path}/scripts/Summary_merge.py",
            f"03-CellReadsSummary/{sample}-CellReads_Summary.xls",
            f"02-STARsolo/{sample}-Solo.out/GeneFull_Ex50pAS/Summary.csv",
            f"03-CellReadsSummary/{sample}-Summary.xls",
            sample
        ])
    print("细胞读数汇总完成")

    # ------------------------ 步骤4：矩阵文件更名 ------------------------
    print("\n===== 开始 矩阵文件更名 =====")
    os.makedirs("04-CB3_Matrix/raw", exist_ok=True)
    os.makedirs("04-CB3_Matrix/filtered", exist_ok=True)
    
    for sample in samples:
        # 复制文件
        subprocess.run([
            "cp", "-r",
            f"02-STARsolo/{sample}-Solo.out/GeneFull_Ex50pAS/raw",
            f"04-CB3_Matrix/raw/{sample}_raw_feature_bc_matrix"
        ])
        subprocess.run([
            "cp", "-r",
            f"02-STARsolo/{sample}-Solo.out/GeneFull_Ex50pAS/filtered",
            f"04-CB3_Matrix/filtered/{sample}_filtered_feature_bc_matrix"
        ])
        
        # 压缩文件
        pigz_cmds = []
        for dir_type in ["raw", "filtered"]:
            for f in ["matrix.mtx", "features.tsv", "barcodes.tsv"]:
                path = f"04-CB3_Matrix/{dir_type}/{sample}_{dir_type}_feature_bc_matrix/{f}"
                pigz_cmds.append(subprocess.Popen(["pigz", "-p", "1", path]))
        
        for cmd in pigz_cmds:
            cmd.wait()
    print("矩阵文件处理完成")

    # ------------------------ 步骤5：合并子文库汇总 ------------------------
    print("\n===== 开始 子文库汇总合并 =====")
    os.makedirs("05-Combined", exist_ok=True)
    summary_files = " ".join(glob.glob("03-CellReadsSummary/*-Summary.xls"))
    subprocess.run(
        f"python {args.script_path}/scripts/Combined_Sample.py "
        f"{summary_files} 05-Combined/1-Results_Summary.xls",
        shell=True,
        check=True
    )
    print("子文库汇总完成")

    # ------------------------ 步骤6：合并CB3子文库 ------------------------
    print("\n===== 开始 CB3子文库合并 =====")
    cellreads_files = " ".join(glob.glob("02-STARsolo/*-Solo.out/GeneFull_Ex50pAS/CellReads.stats"))
    matrix_dirs = " ".join(glob.glob("04-CB3_Matrix/filtered/*_filtered_feature_bc_matrix"))
    
    subprocess.run(
        f"python {args.script_path}/scripts/Cell3_merge.py "
        f"--SplitSample {args.SplitLibrary} "
        f"--inputfileSummary 05-Combined/1-Results_Summary.xls "
        f"--inputfileCellReads {cellreads_files} "
        f"--inputfileMatrixDir {matrix_dirs} "
        f"--outputfile1 05-Combined/2-Cell_Summary.xls "
        f"--outputfile2 05-Combined/2-Cell3_Merge_Summary.xls ",
        shell=True,
        check=True
    )
    print("CB3合并完成")

    # ------------------------ 步骤7：CB2混样拆分 ------------------------
    print("\n===== 开始 CB2混样拆分 =====")
    os.makedirs("06-CB2_Sample", exist_ok=True)
    matrix_paths = " ".join(glob.glob("04-CB3_Matrix/filtered/*_filtered_feature_bc_matrix"))
    
    subprocess.run(
        f"python {args.script_path}/scripts/CellNum.py "
        f"--inputfileCellMatrix {matrix_paths} "
        f"--inputfileBarcodelist {args.script_path}/Barcode/CB2.list "
        f"--splitCB {args.splitCB} "
        f"--splitSample {args.splitSample} "
        f"--outCellNum 4-Sample_CellNum.tsv "
        f"--outCellList {args.splitSample_EstimatedCell_list} ",
        shell=True,
        check=True
    )
    print("CB2拆分完成")

    # ------------------------ 步骤8：生成H5ad文件 ------------------------
    print("\n===== 开始 生成H5ad文件 =====")
    os.makedirs("04-CB3_Matrix/filtered_AddName", exist_ok=True)
    os.makedirs("04-CB3_Matrix/raw_AddName", exist_ok=True)
    
    # 修改Barcode名称
    subprocess.run(
        f"python {args.script_path}/scripts/ModifyBarcodeName.py "
        f"--ModifyBarcodeName True "
        f"--inputMatrixDir 04-CB3_Matrix/filtered "
        f"--outputMatrixDir 04-CB3_Matrix/filtered_AddName ",
        shell=True,
        check=True
    )
    subprocess.run(
        f"python {args.script_path}/scripts/ModifyBarcodeName.py "
        f"--ModifyBarcodeName True "
        f"--inputMatrixDir 04-CB3_Matrix/raw "
        f"--outputMatrixDir 04-CB3_Matrix/raw_AddName ",
        shell=True,
        check=True
    )

    # 生成H5ad
    subprocess.run(
        f"python {args.script_path}/scripts/Scanpy_T.py "
        f"--inputMatrixDir 04-CB3_Matrix/filtered_AddName/ "
        f"--inputfileBarcodelist {args.script_path}/Barcode/CB2.list "
        f"--splitCB {args.splitCB} "
        f"--splitSample {args.splitSample} ",
        shell=True,
        check=True
    )
    print("H5ad生成完成")

    # ------------------------ 步骤9：转换为10X矩阵 ------------------------
    print("\n===== 开始 转换为10X矩阵 =====")
    os.makedirs("07-CB2_Matrix", exist_ok=True)
    samples_split = args.splitSample.split()
    
    processes = []
    for sample in samples_split:
        p = subprocess.Popen([
            "python", f"{args.script_path}/scripts/H5adto10X.py",
            sample,
            f"07-CB2_Matrix/{sample}_filtered_feature_bc_matrix"
        ])
        processes.append(p)
    
    for p in processes:
        p.wait()
    print("10X矩阵转换完成")

    # ------------------------ 步骤10：添加SampleID到汇总表 ------------------------
    print("\n===== 开始 添加SampleID到汇总表 =====")
    # 预处理表格
    file_path = "05-Combined/2-Cell_Summary.xls"
    if not os.path.exists(file_path):
        print(f"文件 {file_path} 不存在，跳过操作")
        sys.exit(0)

    # 读取文件第一行并处理
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            first_line = f.readline().strip()  # 读取第一行并去除首尾空白
    except Exception as e:
        print(f"读取文件失败: {str(e)}")
        sys.exit(1)
    
    # 检查第一行第一列是否为"CB"
    if first_line:
        first_column = first_line.split('\t')[0]  # 按制表符分割后取第一列
        if first_column == "CB":
            print("第一列已为'CB'，跳过插入表头操作")
            sys.exit(0)

    # 如果不满足条件，执行sed插入表头命令
    try:
        subprocess.run(
            [
                "sed", "-i",
                "1iCB\tcbMatch\tcbPerfect\tcbMMunique\tcbMMmultiple\tgenomeU\tgenomeM\tfeatureU\tfeatureM\texonic\tintronic\texonicAS\tintronicAS\tmito\tcountedU\tcountedM\tnUMIunique\tnGenesUnique\tnUMImulti\tnGenesMulti",
                file_path
            ],
            check=True,
            shell=False  # 推荐使用shell=False确保安全
        )
        print(f"已向 {file_path} 插入表头行")
    except subprocess.CalledProcessError as e:
        print(f"执行sed命令失败: {str(e)}")
        sys.exit(1)
    # subprocess.run([
    #     "sed", "-i",
    #     "1iCB\tcbMatch\tcbPerfect\tcbMMunique\tcbMMmultiple\tgenomeU\tgenomeM\tfeatureU\tfeatureM\texonic\tintronic\texonicAS\tintronicAS\tmito\tcountedU\tcountedM\tnUMIunique\tnGenesUnique\tnUMImulti\tnGenesMulti",
    #     "05-Combined/2-Cell_Summary.xls"
    # ])
    subprocess.run(
        f"cut -f 1,2,6,7,15,17,18 "
        f"05-Combined/2-Cell_Summary.xls "
        f"> 05-Combined/2-Cell_Summary_Modify.xls",
        shell=True,
        check=True
    )
    
    # 运行AddSampleID.py
    subprocess.run(
        f"python {args.script_path}/scripts/AddSampleID.py "
        f"--inputfile 05-Combined/2-Cell_Summary_Modify.xls "
        f"--barcodelist {args.script_path}/Barcode/CB2.list "
        f"--splitCB {args.splitCB} "
        f"--splitSample {args.splitSample} "
        f"--output 05-Combined/2-Cell_Summary_withSampleID.xls ",
        shell=True,
        check=True
    )
    print("SampleID添加完成")

    # ------------------------ 步骤11：CB2拆分排序 ------------------------
    print("\n===== 开始 CB2拆分排序 =====")
    sample_ids = args.SampleName.split()
    for sample_id in sample_ids:
        cellreads_files = " ".join(glob.glob(f"02-STARsolo/{sample_id}*-Solo.out/GeneFull_Ex50pAS/CellReads.stats"))
        
        # 运行AddCB2_ID.py
        subprocess.run(
            f"python {args.script_path}/scripts/AddCB2_ID.py "
            f"--inputfile {cellreads_files} "
            f"--outputfile1 05-Combined/2-Cell_AddCB2ID_{sample_id}_Summary.xls "
            f"--group_col CB2_ID "
            f"--stat_cols cbMatch genomeU countedU nUMIunique nGenesUnique "
            f"--outputfile2 05-Combined/CB2ID_{sample_id}_Summary.xls ",
            shell=True,
            check=True
        )
        
        # 运行CB2_sort.py
        subprocess.run(
            f"python {args.script_path}/scripts/CB2_sort.py "
            f"05-Combined/CB2ID_{sample_id}_Summary.xls "
            f"05-Combined/CB2ID_{sample_id}_Summary_sort.xls ",
        shell=True,
        check=True
        )
    print("CB2排序完成")

    # ------------------------ 步骤12：混样统计汇总 ------------------------
    print("\n===== 开始 混样统计汇总 =====")
    os.makedirs("08-Sample_Summary", exist_ok=True)
    cellreads_all = " ".join(glob.glob("02-STARsolo/*-Solo.out/GeneFull_Ex50pAS/CellReads.stats"))
    
    subprocess.run(
        f"python {args.script_path}/scripts/SampleID_split.py "
        f"--inputfileCellReads {cellreads_all} "
        f"--inputfileBarcodelist {args.script_path}/Barcode/CB2.list "
        f"--inputfileCellSummaryFilter 05-Combined/2-Cell_Summary_withSampleID.xls "
        f"--splitCB {args.splitCB} "
        f"--splitSample {args.splitSample} "
        f"--inputfileCellMatrix {args.splitSample_EstimatedCellMatrix} "
        f"--outputfileCellSummaryRaw 08-Sample_Summary/2-Cell_Summary_withSampleID_AllCB.xls "
        f"--outputfileSampleSummary 08-Sample_Summary/5-Sample_Summary.xls ",
        shell=True,
        check=True
    )
    print("混样统计完成")

    # ------------------------ 步骤13：生成报告 ------------------------
    subprocess.run(f"python {args.script_path}/scripts/Report_html.py", shell=True, check=True)
    print("\n===== 10K 超高通量单细胞转录组分析报告已生成！ =====")
    print("\n===== 10K 超高通量单细胞转录组分析全部流程执行完毕！ =====")
if __name__ == "__main__":
    main()