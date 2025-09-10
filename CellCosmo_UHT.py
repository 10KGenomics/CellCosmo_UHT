import os
import glob
import subprocess
import argparse
import sys
from argparse import RawTextHelpFormatter
import shlex

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="10K Ultra-High-Throughput Single-Cell Transcriptome Analysis Pipeline",
        formatter_class=RawTextHelpFormatter
    )
    # Basic parameters
    parser.add_argument("--script_path", required=True, help="Path to pipeline scripts, e.g., /path/CellCosmo_UHT")
    parser.add_argument("--SampleName", required=True, help="Library name(s), space-separated for multiple samples, e.g., 'WT1-1 WT1-2 WT1-3'")
    parser.add_argument("--Rawdata_path", required=True, help="Directory containing Sample_R1.fastq.gz and Sample_R2.fastq.gz files, e.g., './0-data'")
    
    # Fastq processing parameters
    parser.add_argument("--CB3_Num", required=True, help="CB3 index selection range, e.g., 1-20")
    
    # STARsolo parameters
    parser.add_argument("--STARindex", required=True, help="Reference genome STAR index path, e.g., /path/GRCh38_index")
    parser.add_argument("--TopCells", type=int, required=True, help="Forced cell count per CB3 sub-library based on experimental design, e.g., 2000")
    parser.add_argument("--STARsoloThreads", type=int, required=True, help="Thread count per STARsolo job, e.g., 4")
    parser.add_argument("--SplitcodeNum", type=int, required=True, help="Number of splitcode libraries, e.g., 4")
    parser.add_argument("--nFastqs", type=int, required=True, help="Number of fastq files provided, 2 for paired-end, 1 for single-end, e.g., 2")
    parser.add_argument("--soloFeatures", default='GeneFull GeneFull_Ex50pAS Velocyto', 
                        help="Feature types to extract from single-cell data, e.g., 'GeneFull GeneFull_Ex50pAS Velocyto'")
    parser.add_argument("--soloBAM", choices=['True', 'False'], default='True',
                        help="Whether to generate BAM files; default: True(--outSAMtype BAM SortedByCoordinate --outSAMattributes CR UR)")
    parser.add_argument("--soloUMIlen", type=int, default=8, help="UMI length (default: 8)")
    parser.add_argument("--STARsolo_param", help="Additional STARsolo parameters, e.g., '--outReadsUnmapped Fastx'")
    parser.add_argument("--Summary", default='summary.txt', help="Output summary file (default: 'summary.txt')")
    parser.add_argument("--Mapping", default='mapping.txt', help="Mapping output file (default: 'mapping.txt')")
    parser.add_argument("--LogFile", default='log.txt', help="Log file (default: 'run_log.txt')")
    
    # Sample demultiplexing parameters
    parser.add_argument("--outRaw", choices=['True', 'False'], default='False', help="Whether to output Raw matrices & H5ad files; default: False")
    parser.add_argument("--SplitLibrary", default='False', 
                        help="For multiple fastq.gz files from a single large sample. Set to False for single library, or provide batch name for multiple libraries; default: False")
    parser.add_argument("--splitCB", required=True, 
                        help="Demultiplexing classification rules: x-y (1-64 65-128 129-192) or [start,end,step] ([1,89,8] [2,90,8] [3,91,8]) or x,y,z w,a,b (1,2,3,4 5,6,7,8 9,10,11,12) (space-separated for different samples, comma-separated within samples)")
    
    # Add splitSample parameter and related parameters
    parser.add_argument("--splitSample", required=True, 
                        help="List of sample names after demultiplexing, e.g., 'A B C'")
    parser.add_argument("--splitSample_EstimatedCell_list", 
                        help="[Optional] List file of estimated cell counts for each sample")
    parser.add_argument("--splitSample_RawCell_list", 
                        help="[Optional] List file of raw cell counts for each sample")
    parser.add_argument("--splitSample_EstimatedCellMatrix", 
                        help="[Optional] Matrix file prefix for each sample")
    
    args = parser.parse_args()
    
    # Automatically generate values for related parameters
    samples = args.splitSample.split()
    if not args.splitSample_EstimatedCell_list:
        args.splitSample_EstimatedCell_list = " ".join(f"4-EstimatedCell_{s}.list" for s in samples)
    
    if not args.splitSample_RawCell_list:
        args.splitSample_RawCell_list = " ".join(f"4-RawCell_{s}.list" for s in samples)
    
    if not args.splitSample_EstimatedCellMatrix:
        args.splitSample_EstimatedCellMatrix = " ".join(f"{s}_filtered_feature_bc_matrix" for s in samples)
    
    # Print generated parameters for verification
    print(f"Generated parameter values:")
    print(f"splitSample_EstimatedCell_list: {args.splitSample_EstimatedCell_list}")
    print(f"splitSample_RawCell_list: {args.splitSample_RawCell_list}")
    print(f"splitSample_EstimatedCellMatrix: {args.splitSample_EstimatedCellMatrix}")
    
    # Initialize environment
    subprocess.run(f"mkdir -p log/", shell=True)
    os.makedirs("log", exist_ok=True)
    print("===== Environment initialization completed =====")

    print("\n===== Starting Generate lib processing =====")
    fastq_files = glob.glob("log/10K_lib")
    if not fastq_files:
        cmd = [
            f"nohup bash {args.script_path}/scripts/Generate_lib.sh",
            f"--CB3_Num '{args.CB3_Num}'",
            f"--outCB3_lib 10K_lib"
        ]
        subprocess.run(" ".join(cmd), shell=True)
        print("Generate lib processing initiated")
    else:
        print("Existing Generate lib results detected, skipping this step")

    print("\n===== Starting STARsolo analysis =====")
    star_logs = glob.glob("02-STARsolo/star_outs_lib*/Log.final.out")
    if not star_logs:
        # Process input files for multiple samples
        sample_names = args.SampleName.split()  # Split sample names
        
        # Build input file list in R2 R1 R2 R1 order
        input_files = []
        for sample in sample_names:
            r2_files = glob.glob(f"{args.Rawdata_path}/{sample}_*2.f*q.gz")
            r1_files = glob.glob(f"{args.Rawdata_path}/{sample}_*1.f*q.gz")
            
            if not r2_files or not r1_files:
                print(f"Error: Cannot find fastq files for sample {sample}")
                sys.exit(1)
                
            # Add to list in R2 R1 order
            input_files.extend([r2_files[0], r1_files[0]])
        
        # Convert file list to string
        input_files_str = " ".join(input_files)

        if args.soloBAM == 'True':
            bam_option = '--outSAMtype BAM SortedByCoordinate --outSAMattributes CR UR'
        else:
            bam_option = '--outSAMtype None'
        
        cmd_base = [
            "nohup", "bash", f"{args.script_path}/scripts/STARsolo.sh",
            "--STARindex", args.STARindex,
            "--CB3_Num", args.CB3_Num,
            "--STARsoloThreads", str(args.STARsoloThreads),
            "--TopCell", str(args.TopCells),
            "--SplitcodeNum", str(args.SplitcodeNum),
            "--nFastqs", str(args.nFastqs),
            "--config", f"{args.script_path}/scripts/config_Demultiplexing.txt",
            "--inputFiles", input_files_str,
            "--soloFeatures", args.soloFeatures,
            "--soloBAM_param", bam_option,
            "--soloUMIlen", str(args.soloUMIlen),
            "--SplitcodeLib", "log/10K_lib",
            "--Summary", args.Summary,
            "--Mapping", args.Mapping,
            "--LogFile", args.LogFile
        ]

        if args.STARsolo_param:
            cmd_base.extend(["--STARsolo_param", args.STARsolo_param])
    
        print(f"STARsolo initiated, see log: {args.LogFile}")
        with open(args.LogFile, "w") as log_file:
            subprocess.run(cmd_base, stdout=log_file, stderr=subprocess.STDOUT)
    else:
        print("Existing STARsolo results detected, skipping this step")

    print("\n===== Starting cell reads summary =====")
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
    print(f"Found {len(samples)} samples: {samples}")

    for sample in samples:
        print(f"\nProcessing sample: {sample}")
    
        cellreads_path = f"02-STARsolo/star_outs_lib{sample}/Solo.out/GeneFull_Ex50pAS/CellReads.stats"
        summary_path = f"02-STARsolo/star_outs_lib{sample}/Solo.out/GeneFull_Ex50pAS/Summary.csv"
    
        if not os.path.exists(cellreads_path):
            print(f"Warning: File {cellreads_path} does not exist, skipping")
            continue  
        if not os.path.exists(summary_path):
            print(f"Warning: File {summary_path} does not exist, skipping")
            continue
    
        print(f"Processing: {cellreads_path}")
        try:
            subprocess.run([
                "python", f"{args.script_path}/scripts/CellReads_Summary.py",
                cellreads_path,
                f"03-CellReadsSummary/{sample}-CellReads_Summary.xls",
                sample
            ], check=True)
        except subprocess.CalledProcessError as e:
            print(f"CellReads_Summary.py execution failed: {str(e)}")
    
        print(f"Processing: {summary_path}")
        try:
            subprocess.run([
                "python", f"{args.script_path}/scripts/Summary_merge.py",
                f"03-CellReadsSummary/{sample}-CellReads_Summary.xls",
                summary_path,
                f"03-CellReadsSummary/{sample}-Summary.xls",
                sample
            ], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Summary_merge.py execution failed: {str(e)}")

    print("Cell reads summary completed")

    print("\n===== Starting matrix file renaming =====")
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
    print("Matrix file processing completed")

    print("\n===== Starting sub-library summary and merging =====")
    os.makedirs("05-Combined", exist_ok=True)
    summary_files = " ".join(glob.glob("03-CellReadsSummary/*-Summary.xls"))
    subprocess.run(
        f"python {args.script_path}/scripts/Combined_Sample.py "
        f"{summary_files} 05-Combined/1-Results_Summary.xls",
        shell=True,
        check=True
    )
    print("Sub-library summary completed")

    print("\n===== Starting sub-library merging =====")
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
    print("Merging completed")

    print("\n===== Starting sample demultiplexing =====")
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
    print("Filtered demultiplexing completed")

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
        print("Raw demultiplexing completed")

    print("\n===== Starting H5ad file generation =====")
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
    print("Filtered H5ad generation completed")

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
        print("Raw H5ad generation completed")

    print("\n===== Starting 10X matrix conversion =====")
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
    print("10X filtered matrix conversion completed")

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
        print("10X raw matrix conversion completed")

    print("\n===== Starting SampleID addition to summary table =====")
    file_path = "05-Combined/2-Cell_Summary.xls"
    if not os.path.exists(file_path):
        print(f"ERROR: File {file_path} does not exist, cannot proceed", file=sys.stderr)
        sys.exit(1)

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            first_line = f.readline().strip()  
    except Exception as e:
        print(f"ERROR: File reading failed: {str(e)}", file=sys.stderr)
        sys.exit(1)
    
    insert_header = True
    if first_line:
        first_column = first_line.split('\t')[0] 
        if first_column == "CB":
            print("First column is already 'CB', skipping header insertion")
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
            print(f"Header row inserted into {file_path}")
        except subprocess.CalledProcessError as e:
            print(f"ERROR: sed command execution failed: {str(e)}", file=sys.stderr)
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
        print("Generated 05-Combined/2-Cell_Summary_Modify.xls")
    except subprocess.CalledProcessError as e:
        print(f"ERROR: cut command execution failed: {str(e)}", file=sys.stderr)
        sys.exit(1)

    modify_file = "05-Combined/2-Cell_Summary_Modify.xls"
    if not os.path.exists(modify_file) or os.path.getsize(modify_file) == 0:
        print(f"ERROR: Generated {modify_file} does not exist or is empty", file=sys.stderr)
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
        
        print(f"Executing command: {' '.join(cmd)}")
        
        subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        print("SampleID addition completed")
    except subprocess.CalledProcessError as e:
        print(f"ERROR: AddSampleID.py execution failed: {e.stderr}", file=sys.stderr)
        sys.exit(1)

    result_file = "05-Combined/2-Cell_Summary_withSampleID.xls"
    if not os.path.exists(result_file) or os.path.getsize(result_file) == 0:
        print(f"ERROR: Result file {result_file} was not generated or is empty", file=sys.stderr)
        sys.exit(1)

    print("\n===== Starting demultiplexing statistics summary =====")
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
        print(f"WARNING: Potential issue with demultiplexing statistics summary step: {str(e)}")

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

    print("Demultiplexing statistics completed")

    try:
        subprocess.run(f"python {args.script_path}/scripts/Report_html.py", shell=True, check=True)
        print("\n===== 10K Ultra-High-Throughput Single-Cell Transcriptome Analysis Report Generated! =====")
    except Exception as e:
        print(f"WARNING: Potential issue with report generation step: {str(e)}")

    print("\n===== 10K Ultra-High-Throughput Single-Cell Transcriptome Analysis Pipeline Execution Completed! =====")

if __name__ == "__main__":
    main()