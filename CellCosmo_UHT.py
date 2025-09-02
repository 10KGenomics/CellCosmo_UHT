import os
import glob
import subprocess
import argparse
import sys
from argparse import RawTextHelpFormatter
import shlex

def main():
    # 解析命令行参数
    parser = argparse.ArgumentParser(
        description="10K 超高通量单细胞转录组分析流程Pipeline ",
        formatter_class=RawTextHelpFormatter
    )
    # 基础参数
    parser.add_argument("--script_path", required=True, help="pipeline每步脚本存放路径，比如：/path/CellCosmo_UHT")
    parser.add_argument("--SampleName", required=True, help="下机文库名，比如：'WT1-1'")
    parser.add_argument("--Rawdata_path", required=True, help="Sample_R1.fastq.gz、Sample_R2.fastq.gz存放目录，比如：'./0-data'")
    
    # Fastq处理参数
    parser.add_argument("--CB3_Num", required=True, help="CB3选取范围index，比如：1-20")
    
    # STARsolo参数
    parser.add_argument("--STARindex", required=True, help="参考基因组STAR index 路径，比如：/path/GRCh38_index")
    parser.add_argument("--TopCells", type=int, required=True, help="CB3单个子文库的强制细胞数,根据实验信息填写，比如：2000")
    parser.add_argument("--STARsoloThreads", type=int, required=True, help="STARsolo单个任务线程数，比如：4")
    parser.add_argument("--SplitcodeNum", type=int, required=True, help="splitcode lib数，比如：4")
    parser.add_argument("--nFastqs", type=int, required=True, help="提供fastq文件数，双端测序=2，单端=1，比如：2")
    parser.add_argument("--soloFeatures", default='GeneFull GeneFull_Ex50pAS Velocyto', 
                        help="需要从单细胞数据中提取的特征类型,比如：'GeneFull GeneFull_Ex50pAS Velocyto'")
    parser.add_argument("--soloUMIlen", type=int, default=8, help="UMI长度 (默认: 8)")
    parser.add_argument("--STARsolo_param", help="STARsolo增加参数，比如：'--outReadsUnmapped Fastx'")
    parser.add_argument("--Summary", default='summary.txt', help="要输出文件 (默认: 'summary.txt')")
    parser.add_argument("--Mapping", default='mapping.txt', help="映射输出文件 (默认: 'mapping.txt')")
    parser.add_argument("--LogFile", default='log.txt', help="日志文件 (默认: 'run_log.txt')")
    
    # 混样拆分参数
    parser.add_argument("--outRaw", choices=['True', 'False'], default='False', help="是否输出Raw 矩阵&H5ad；默认：False")
    parser.add_argument("--SplitLibrary", default='False', 
                        help="下机多组fastq.gz,实际1个大样本。若单个下机文库名，则填写False；若多个下机文库名，则填写批次名；默认：False")
    parser.add_argument("--splitCB", required=True, 
                        help="混样CB2分类规则：x-y（1-64 65-128 129-192）或[起始值,终止值,公差值]（[1,89,8] [2,90,8] [3,91,8]）或x,y,z w,a,b（1,2,3,4 5,6,7,8 9,10,11,12） （以空格分隔不同样本，逗号分隔样本内部CB2的元素）")
    
    # 添加splitSample参数及其关联参数
    parser.add_argument("--splitSample", required=True, 
                        help="混样拆分后的样本名列表,比如：'A B C'")
    parser.add_argument("--splitSample_EstimatedCell_list", 
                        help="[可选] 各样本预估细胞数列表文件")
    parser.add_argument("--splitSample_RawCell_list", 
                        help="[可选] 各样本原始细胞数列表文件")
    parser.add_argument("--splitSample_EstimatedCellMatrix", 
                        help="[可选] 各样本矩阵文件前缀")
    
    args = parser.parse_args()
    
    # 自动生成关联参数的值
    samples = args.splitSample.split()
    if not args.splitSample_EstimatedCell_list:
        args.splitSample_EstimatedCell_list = " ".join(f"4-EstimatedCell_{s}.list" for s in samples)
    
    if not args.splitSample_RawCell_list:
        args.splitSample_RawCell_list = " ".join(f"4-RawCell_{s}.list" for s in samples)
    
    if not args.splitSample_EstimatedCellMatrix:
        args.splitSample_EstimatedCellMatrix = " ".join(f"{s}_filtered_feature_bc_matrix" for s in samples)
    
    print(f"生成的参数值:")
    print(f"splitSample_EstimatedCell_list: {args.splitSample_EstimatedCell_list}")
    print(f"splitSample_RawCell_list: {args.splitSample_RawCell_list}")
    print(f"splitSample_EstimatedCellMatrix: {args.splitSample_EstimatedCellMatrix}")
    
    # parser.add_argument("--CB2_list", required=True, help="CB2 barcode列表文件绝对路径，比如：/path/CellCosmo_UHT/Barcode/CB2.list ")
    # parser.add_argument("--condaName", required=True, help="conda环境名称，比如：CellCosmo_UHT")
    # parser.add_argument("--config", required=True, help="splitcode配置文件，比如：/path/config_Demultiplexing.txt")
    # parser.add_argument("--splitSample_EstimatedCell_list", required=True, help="各样本预估细胞数列表文件,与splitSample对应数目，比如：'4-EstimatedCell_A.list 4-EstimatedCell_B.list  4-EstimatedCell_C.list'")
    # parser.add_argument("--splitSample_RawCell_list", required=True, help="各样本原始细胞数列表文件,与splitSample对应数目，比如：'4-RawCell_A.list 4-RawCell_B.list  4-RawCell_C.list'")
    # parser.add_argument("--splitSample_EstimatedCellMatrix", required=True, help="各样本矩阵文件前缀,与splitSample对应数目，比如：'A_filtered_feature_bc_matrix B_filtered_feature_bc_matrix C_filtered_feature_bc_matrix'")

    subprocess.run(f"mkdir -p log/", shell=True)
    os.makedirs("log", exist_ok=True)
    print("===== 环境初始化完成 =====")

    print("\n===== 开始 Generate lib 处理 =====")
    fastq_files = glob.glob("log/10K_lib")
    if not fastq_files:
        cmd = [
            f"nohup bash {args.script_path}/scripts/Generate_lib.sh",
            f"--CB3_Num '{args.CB3_Num}'",
            f"--outCB3_lib 10K_lib"
        ]
        subprocess.run(" ".join(cmd), shell=True)
        print("Generate lib处理已启动 ")
    else:
        print("检测到已存在Generate lib处理结果，跳过此步骤")

    print("\n===== 开始 STARsolo 分析 =====")
    star_logs = glob.glob("02-STARsolo/star_outs_lib*/Log.final.out")
    if not star_logs:
        input_r1 = ' '.join(glob.glob(f"{args.Rawdata_path}/{args.SampleName}_*1.f*q.gz"))
        input_r2 = ' '.join(glob.glob(f"{args.Rawdata_path}/{args.SampleName}_*2.f*q.gz"))

        cmd_base = [
            "nohup", "bash", f"{args.script_path}/scripts/STARsolo.sh",
            "--STARindex", args.STARindex,
            "--CB3_Num", args.CB3_Num,
            "--STARsoloThreads", str(args.STARsoloThreads),
            "--TopCell", str(args.TopCells),
            "--SplitcodeNum", str(args.SplitcodeNum),
            "--nFastqs", str(args.nFastqs),
            "--config", f"{args.script_path}/scripts/config_Demultiplexing.txt",
            "--inputR1", input_r1,
            "--inputR2", input_r2,
            "--soloFeatures", args.soloFeatures,
            "--soloUMIlen", str(args.soloUMIlen),
            "--SplitcodeLib", "log/10K_lib",
            "--Summary", args.Summary,
            "--Mapping", args.Mapping,
            "--LogFile", args.LogFile
        ]

        if args.STARsolo_param:
            cmd_base.extend(["--STARsolo_param", args.STARsolo_param])
    
        print(f"STARsolo已启动，日志见: {args.LogFile}")
        with open(args.LogFile, "w") as log_file:
            subprocess.run(cmd_base, stdout=log_file, stderr=subprocess.STDOUT)
    else:
        print("检测到已存在STARsolo结果，跳过此步骤")

    print("\n===== 开始 细胞读数汇总 =====")
    os.makedirs("03-CellReadsSummary", exist_ok=True)

    star_dirs = glob.glob("02-STARsolo/star_outs_lib*")
    samples = []
    for d in star_dirs:
        if "star_outs_lib" in d:
            parts = d.split("star_outs_lib")
            if len(parts) > 1:
                sample_id = parts[-1].split("/")[0]
                if sample_id.isdigit():
                    samples.append(sample_id)

    samples.sort(key=int)
    print(f"找到 {len(samples)} 个样本: {samples}")

    for sample in samples:
        print(f"\n处理样本: {sample}")
    
        cellreads_path = f"02-STARsolo/star_outs_lib{sample}/Solo.out/GeneFull_Ex50pAS/CellReads.stats"
        summary_path = f"02-STARsolo/star_outs_lib{sample}/Solo.out/GeneFull_Ex50pAS/Summary.csv"
    
        if not os.path.exists(cellreads_path):
            print(f"警告: 文件 {cellreads_path} 不存在，跳过")
            continue  
        if not os.path.exists(summary_path):
            print(f"警告: 文件 {summary_path} 不存在，跳过")
            continue
    
        print(f"处理: {cellreads_path}")
        try:
            # 运行CellReads_Summary.py
            subprocess.run([
                "python", f"{args.script_path}/scripts/CellReads_Summary.py",
                cellreads_path,
                f"03-CellReadsSummary/{sample}-CellReads_Summary.xls",
                sample
            ], check=True)
        except subprocess.CalledProcessError as e:
            print(f"CellReads_Summary.py 执行失败: {str(e)}")
    
        print(f"处理: {summary_path}")
        try:
            # 运行Summary_merge.py
            subprocess.run([
                "python", f"{args.script_path}/scripts/Summary_merge.py",
                f"03-CellReadsSummary/{sample}-CellReads_Summary.xls",
                summary_path,
                f"03-CellReadsSummary/{sample}-Summary.xls",
                sample
            ], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Summary_merge.py 执行失败: {str(e)}")

    print("细胞读数汇总完成")

    print("\n===== 开始 矩阵文件更名 =====")
    os.makedirs("04-CB3_Matrix/raw", exist_ok=True)
    os.makedirs("04-CB3_Matrix/filtered", exist_ok=True)
    
    for sample in samples:
        subprocess.run([
            "cp", "-r",
            f"02-STARsolo/star_outs_lib{sample}/Solo.out/GeneFull_Ex50pAS/raw",
            f"04-CB3_Matrix/raw/{sample}_raw_feature_bc_matrix"
        ])
        subprocess.run([
            "cp", "-r",
            f"02-STARsolo/star_outs_lib{sample}/Solo.out/GeneFull_Ex50pAS/filtered",
            f"04-CB3_Matrix/filtered/{sample}_filtered_feature_bc_matrix"
        ])
        
        pigz_cmds = []
        for dir_type in ["raw", "filtered"]:
            for f in ["matrix.mtx", "features.tsv", "barcodes.tsv"]:
                path = f"04-CB3_Matrix/{dir_type}/{sample}_{dir_type}_feature_bc_matrix/{f}"
                pigz_cmds.append(subprocess.Popen(["pigz", "-p", "1", path]))
        
        for cmd in pigz_cmds:
            cmd.wait()
    print("矩阵文件处理完成")

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

    print("\n===== 开始 CB3子文库合并 =====")
    cellreads_files = " ".join(glob.glob("02-STARsolo/star_outs_lib*/Solo.out/GeneFull_Ex50pAS/CellReads.stats"))
    matrix_dirs = " ".join(glob.glob("04-CB3_Matrix/filtered/*_filtered_feature_bc_matrix"))
    
    subprocess.run(
        f"python {args.script_path}/scripts/Cell3_merge.py "
        f"--SplitSample {args.SplitLibrary} "
        f"--inputfileSummary 05-Combined/1-Results_Summary.xls "
        f"--inputfileCellReads {cellreads_files} "
        f"--inputfileMatrixDir {matrix_dirs} "
        f"--inputfileSplitSummary   {args.Summary} "
        f"--outputfile1 05-Combined/2-Cell_Summary.xls "
        f"--outputfile2 05-Combined/2-Cell3_Merge_Summary.xls ",
        shell=True,
        check=True
    )

    print("CB3合并完成")

    print("\n===== 开始 CB2混样拆分 =====")
    os.makedirs("06-CB2_Sample", exist_ok=True)

    matrix_paths_filtered = " ".join(glob.glob("04-CB3_Matrix/filtered/*_filtered_feature_bc_matrix"))
    subprocess.run(
        f"python {args.script_path}/scripts/CellNum.py "
        f"--inputfileCellMatrix {matrix_paths_filtered} "
        f"--inputfileBarcodelist {args.Mapping} "
        f"--splitCB {args.splitCB} "
        f"--splitSample {args.splitSample} "
        f"--outCellNum 4-Sample_CellNum_filtered.tsv "
        f"--outCellList {args.splitSample_EstimatedCell_list} ",
        shell=True,
        check=True
    )
    print("filtered CB2拆分完成")

    if args.outRaw == 'True':
        matrix_paths_raw = " ".join(glob.glob("04-CB3_Matrix/raw/*_raw_feature_bc_matrix"))
        subprocess.run(
            f"python {args.script_path}/scripts/CellNum.py "
            f"--inputfileCellMatrix {matrix_paths_raw} "
            f"--inputfileBarcodelist {args.Mapping} "
            f"--splitCB {args.splitCB} "
            f"--splitSample {args.splitSample} "
            f"--outCellNum 4-Sample_CellNum_raw.tsv "
            f"--outCellList {args.splitSample_RawCell_list} ",
            shell=True,
            check=True
        )
        print("raw CB2拆分完成")

    print("\n===== 开始 生成H5ad文件 =====")
    os.makedirs("04-CB3_Matrix/filtered_AddName", exist_ok=True)
    os.makedirs("04-CB3_Matrix/raw_AddName", exist_ok=True)
    
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

    subprocess.run(
        f"python {args.script_path}/scripts/Scanpy_T.py "
        f"--inputMatrixDir 04-CB3_Matrix/filtered_AddName "
        f"--matrixType filtered "
        f"--inputfileBarcodelist  {args.Mapping} "
        f"--splitCB {args.splitCB} "
        f"--splitSample {args.splitSample} ",
        shell=True,
        check=True
    )
    print("filtered H5ad生成完成")

    if args.outRaw == 'True':
        subprocess.run(
            f"python {args.script_path}/scripts/Scanpy_T.py "
            f"--inputMatrixDir 04-CB3_Matrix/raw_AddName "
            f"--matrixType raw "
            f"--inputfileBarcodelist  {args.Mapping} "
            f"--splitCB {args.splitCB} "
            f"--splitSample {args.splitSample} ",
            shell=True,
            check=True
        )
        print("Raw H5ad生成完成")

    print("\n===== 开始 转换为10X矩阵 =====")
    os.makedirs("07-CB2_Matrix", exist_ok=True)
    samples_split = args.splitSample.split()  

    processes_filtered = []  
    for sample in samples_split:
        p = subprocess.Popen([
            "python", f"{args.script_path}/scripts/H5adto10X.py",
            f"06-CB2_Sample/filtered/{sample}_filtered.h5ad",
            f"07-CB2_Matrix/{sample}_filtered_feature_bc_matrix"
        ])
        processes_filtered.append(p)
    for p in processes_filtered:
        p.wait()
    print("10X filtered 矩阵转换完成")

    if args.outRaw == 'True':
        processes_raw = []  
        for sample in samples_split:  
            p = subprocess.Popen([
                "python", f"{args.script_path}/scripts/H5adto10X.py",
                f"06-CB2_Sample/raw/{sample}_raw.h5ad",
                f"07-CB2_Matrix/{sample}_raw_feature_bc_matrix"
            ])
            processes_raw.append(p)
        for p in processes_raw:
            p.wait()
        print("10X raw 矩阵转换完成")

    print("\n===== 开始 添加SampleID到汇总表 =====")
    # 预处理表格
    file_path = "05-Combined/2-Cell_Summary.xls"
    if not os.path.exists(file_path):
        print(f"ERROR: 文件 {file_path} 不存在，无法继续", file=sys.stderr)
        sys.exit(1)

    # 读取文件第一行并处理
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            first_line = f.readline().strip()  
    except Exception as e:
        print(f"ERROR: 读取文件失败: {str(e)}", file=sys.stderr)
        sys.exit(1)
    
    insert_header = True
    if first_line:
        first_column = first_line.split('\t')[0]  
        if first_column == "CB":
            print("第一列已为'CB'，跳过插入表头操作")
            insert_header = False

    if insert_header:
        try:
            subprocess.run(
                [
                    "sed", "-i",
                    "1iCB\tcbMatch\tcbPerfect\tcbMMunique\tcbMMmultiple\tgenomeU\tgenomeM\tfeatureU\tfeatureM\texonic\tintronic\texonicAS\tintronicAS\tmito\tcountedU\tcountedM\tnUMIunique\tnGenesUnique\tnUMImulti\tnGenesMulti",
                    file_path
                ],
                check=True,
                shell=False
            )
            print(f"已向 {file_path} 插入表头行")
        except subprocess.CalledProcessError as e:
            print(f"ERROR: 执行sed命令失败: {str(e)}", file=sys.stderr)
            sys.exit(1)

    try:
        subprocess.run(
            f"cut -f 1,2,6,7,15,17,18 "
            f"{file_path} "
            f"> 05-Combined/2-Cell_Summary_Modify.xls",
            shell=True,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        print("已生成 05-Combined/2-Cell_Summary_Modify.xls")
    except subprocess.CalledProcessError as e:
        print(f"ERROR: 执行cut命令失败: {str(e)}", file=sys.stderr)
        sys.exit(1)

    modify_file = "05-Combined/2-Cell_Summary_Modify.xls"
    if not os.path.exists(modify_file) or os.path.getsize(modify_file) == 0:
        print(f"ERROR: 生成的 {modify_file} 不存在或为空", file=sys.stderr)
        sys.exit(1)

    try:
        cmd = [
            "python", f"{args.script_path}/scripts/AddSampleID.py",
            "--inputfile", "05-Combined/2-Cell_Summary_Modify.xls",
            "--barcodelist", args.Mapping,
            "--splitCB", args.splitCB,
            "--splitSample", args.splitSample,
            "--output", "05-Combined/2-Cell_Summary_withSampleID.xls"
        ]
        
        print(f"执行命令: {' '.join(cmd)}")
        
        subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        print("SampleID添加完成")
    except subprocess.CalledProcessError as e:
        print(f"ERROR: 执行AddSampleID.py失败: {e.stderr}", file=sys.stderr)
        sys.exit(1)

    result_file = "05-Combined/2-Cell_Summary_withSampleID.xls"
    if not os.path.exists(result_file) or os.path.getsize(result_file) == 0:
        print(f"ERROR: 结果文件 {result_file} 未生成或为空", file=sys.stderr)
        sys.exit(1)

    print("\n===== 开始 混样统计汇总 =====")
    os.makedirs("08-Sample_Summary", exist_ok=True)
    cellreads_all = " ".join(glob.glob("02-STARsolo/star_outs_lib*/Solo.out/GeneFull_Ex50pAS/CellReads.stats"))
    
    try:
        subprocess.run(
            f"python {args.script_path}/scripts/SampleID_split.py "
            f"--inputfileCellReads {cellreads_all} "
            f"--inputfileBarcodelist {args.Mapping} "
            f"--inputfileCellSummaryFilter 05-Combined/2-Cell_Summary_withSampleID.xls "
            f"--splitCB {args.splitCB} "
            f"--splitSample {args.splitSample} "
            f"--inputfileCellMatrix {args.splitSample_EstimatedCellMatrix if 'splitSample_EstimatedCellMatrix' in args else ''} "
            f"--outputfileCellSummaryRaw 08-Sample_Summary/2-Cell_Summary_withSampleID_AllCB.xls "
            f"--outputfileSampleSummary 08-Sample_Summary/5-Sample_Summary.xls ",
            shell=True,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
    except Exception as e:
        print(f"WARNING: 混样统计汇总步骤可能存在问题: {str(e)}")

    subprocess.run([
        "cp", "-r",
        f"05-Combined/1-Results_Summary.xls",
        f"PerLibrary_Summary.xls"
    ], check=True)

    subprocess.run([
        "cp", "-r",
        f"05-Combined/2-Cell3_Merge_Summary.xls",
        f"MergeLibrary_Summary.xls"
    ], check=True)

    subprocess.run([
        "cp", "-r",
        f"08-Sample_Summary/5-Sample_Summary.xls",
        f"Sample_Summary.xls"
    ], check=True)

    print("混样统计完成")

    try:
        subprocess.run(f"python {args.script_path}/scripts/Report_html.py", shell=True, check=True)
        print("\n===== 10K 超高通量单细胞转录组分析报告已生成！ =====")
    except Exception as e:
        print(f"WARNING: 报告生成步骤可能存在问题: {str(e)}")

    print("\n===== 10K 超高通量单细胞转录组分析全部流程执行完毕！ =====")

if __name__ == "__main__":
    main()
