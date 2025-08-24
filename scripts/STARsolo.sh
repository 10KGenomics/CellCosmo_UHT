#!/bin/bash

# 解析命令行参数
usage() {
    echo "用法: $0 [选项]"
    echo "选项:"
    echo "  --STARindex <路径>          : STAR索引目录 (必需)"
    echo "  --CB3_Num <字符串>          : 样本编号范围 (如 '1-20 24 28') (必需)"
    echo "  --STARsoloThreads <数字>    : STAR线程数 (必需)"
    echo "  --TopCell <数字>            : TopCells过滤值 (必需)"
    echo "  --nFastqs <数字>            : 提供fastq文件数 (必需)"
    echo "  --SplitcodeNum <数字>       : splitcode线程数 (必需)"
    echo "  --config <路径>             : splitcode配置文件 (必需)"
    echo "  --inputR1 <路径>            : 输入R1文件 (必需)"
    echo "  --inputR2 <路径>            : 输入R2文件 (必需)"
    echo "  --soloFeatures <字符串>     : soloFeatures参数 (默认: 'GeneFull GeneFull_Ex50pAS Velocyto')"
    echo "  --soloUMIlen <数字>         : UMI长度 (默认: 8)"
    echo "  --SplitcodeLib <字符串>     : splitcode保留组名 (默认: '10K_lib')"
    echo "  --Summary <文件>            : 摘要输出文件 (默认: 'summary.txt')"
    echo "  --Mapping <文件>            : 映射输出文件 (默认: 'mapping.txt')"
    echo "  --LogFile <文件>            : 日志文件 (默认: 'run_log.txt')"
    echo "  --STARsolo_param <字符串>   : STAR额外参数 (如: '--outReadsUnmapped Fastx')"
    exit 1
}

# 设置默认值
soloFeatures_default="GeneFull GeneFull_Ex50pAS Velocyto"
soloUMIlen_default=8
SplitcodeLib_default="10K_lib"
Summary_default="summary.txt"
Mapping_default="mapping.txt"
LogFile_default="run_log.txt"

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case "$1" in
        --STARindex)
            STARindex="$2"
            shift 2
            ;;
        --CB3_Num)
            CB3_Num="$2"
            shift 2
            ;;
        --STARsoloThreads)
            STARsoloThreads="$2"
            shift 2
            ;;
        --TopCell)
            TopCell="$2"
            shift 2
            ;;
        --nFastqs)
            nFastqs="$2"
            shift 2
            ;;
        --SplitcodeNum)
            SplitcodeNum="$2"
            shift 2
            ;;
        --config)
            config="$2"
            shift 2
            ;;
        --inputR1)
            inputR1="$2"
            shift 2
            ;;
        --inputR2)
            inputR2="$2"
            shift 2
            ;;
        --soloFeatures)
            soloFeatures="$2"
            shift 2
            ;;
        --soloUMIlen)
            soloUMIlen="$2"
            shift 2
            ;;
        --SplitcodeLib)
            SplitcodeLib="$2"
            shift 2
            ;;
        --Summary)
            Summary="$2"
            shift 2
            ;;
        --Mapping)
            Mapping="$2"
            shift 2
            ;;
        --LogFile)
            LogFile="$2"
            shift 2
            ;;
        --STARsolo_param)
            STARsolo_param="$2"
            shift 2
            ;;
        *)
            echo "未知选项: $1"
            usage
            ;;
    esac
done

# 检查必需参数
if [[ -z "$STARindex" || -z "$CB3_Num" || -z "$STARsoloThreads" || -z "$TopCell" || -z "$nFastqs" || 
      -z "$SplitcodeNum" || -z "$config" || -z "$inputR1" || -z "$inputR2" ]]; then
    echo "错误: 缺少必需参数!"
    usage
fi

# 应用默认值（如果未提供）
: "${soloFeatures:=$soloFeatures_default}"
: "${soloUMIlen:=$soloUMIlen_default}"
: "${SplitcodeLib:=$SplitcodeLib_default}"
: "${Summary:=$Summary_default}"
: "${Mapping:=$Mapping_default}"
: "${LogFile:=$LogFile_default}"

# 处理样本范围字符串（将 "1-20 24 28" 转换为数组）
process_ranges() {
    local ranges=($1)
    local result=()
    for range in "${ranges[@]}"; do
        if [[ $range =~ ^([0-9]+)-([0-9]+)$ ]]; then
            start=${BASH_REMATCH[1]}
            end=${BASH_REMATCH[2]}
            for ((i=start; i<=end; i++)); do
                result+=("$i")
            done
        else
            result+=("$range")
        fi
    done
    echo "${result[@]}"
}

# 获取样本数组
samples=($(process_ranges "$CB3_Num"))
# 将数组转换为字符串用于子shell
samples_str="${samples[*]}"

# 主执行块
{
    echo "===$(date '+%F %T')=== 开始运行 ==="
    echo "使用的参数:"
    echo "STAR索引: $STARindex"
    echo "CB3文库index: ${samples[@]}"
    echo "STAR线程: $STARsoloThreads"
    echo "TopCell: $TopCell"
    echo "nFastqs: $nFastqs"
    echo "soloFeatures: $soloFeatures"
    echo "soloUMIlen: $soloUMIlen"
    echo "splitcode线程: $SplitcodeNum"
    echo "配置文件: $config"
    echo "输入文件: R1=$inputR1, R2=$inputR2"
    echo "splitcode保留组: $SplitcodeLib"
    echo "摘要文件: $Summary"
    echo "映射文件: $Mapping"
    echo "额外参数: ${STARsolo_param:-无}"
    echo "日志文件: $LogFile"

    # 使用time命令计时
    /usr/bin/time -v bash -c "
        set -e
        set -o pipefail
        
        # 重建样本数组
        samples_arr=($samples_str)
        echo \"=== 文库数组重建: \${#samples_arr[@]} 个文库 ===\"
        
        # 1. 加载基因组到共享内存
        echo \"正在加载基因组...\"
        if ! STAR --genomeDir \"$STARindex\" --genomeLoad LoadAndExit --outFileNamePrefix log/GenomeLoad; then
            echo \"错误: 基因组加载失败!\"
            exit 1
        fi

        # 2. 创建FIFO管道 (清理旧文件)
        echo \"正在创建FIFO管道...\"
        rm -rf 01-Fastq_Modify_ParaFly/
        rm -rf 02-STARsolo/
        mkdir -p 01-Fastq_Modify_ParaFly/
        mkdir -p 02-STARsolo/
        for s in \"\${samples_arr[@]}\"; do
            # 清理可能存在的旧FIFO
            rm -f 01-Fastq_Modify_ParaFly/\"lib\${s}_R1.fastq\" 01-Fastq_Modify_ParaFly/\"lib\${s}_R2.fastq\" 2>/dev/null
            mkfifo 01-Fastq_Modify_ParaFly/\"lib\${s}_R1.fastq\"
            mkfifo 01-Fastq_Modify_ParaFly/\"lib\${s}_R2.fastq\"
            echo \"创建FIFO: 01-Fastq_Modify_ParaFly/lib\${s}_R1.fastq 和 01-Fastq_Modify_ParaFly/lib\${s}_R2.fastq\"
        done

        # 3. 启动后台STAR进程
        echo \"启动STAR进程... (共\${#samples_arr[@]}个)\"
        for s in \"\${samples_arr[@]}\"; do
            echo \"启动样本\${s}的STAR进程...\"
            STAR \\
                --runThreadN $STARsoloThreads \\
                --genomeDir \"$STARindex\" \\
                --outFileNamePrefix 02-STARsolo/\"star_outs_lib\${s}/\" \\
                --readFilesIn 01-Fastq_Modify_ParaFly/\"lib\${s}_R1.fastq\" 01-Fastq_Modify_ParaFly/\"lib\${s}_R2.fastq\" \\
                --soloUMIlen $soloUMIlen \\
                --soloCBmatchWLtype 1MM \\
                --outSAMattributes CR UR \\
                --soloCBwhitelist None \\
                --soloFeatures $soloFeatures \\
                --outSAMtype BAM SortedByCoordinate \\
                --soloType CB_UMI_Simple \\
                --soloCellReadStats Standard \\
                --soloCellFilter TopCells $TopCell 0.99 10 45000 90000 500 0.01 20000 0.001 10000 \\
                --outFilterScoreMinOverLread 0.33 \\
                --outFilterMatchNminOverLread 0.33 \\
                --genomeLoad LoadAndKeep \\
                --limitBAMsortRAM 1001013131 \\
                $STARsolo_param & 
            echo \"样本\${s}的STAR进程已启动 (PID: \$!)\"
        done

        # 4. 运行splitcode
        echo \"运行splitcode...\"
        splitcode -c \"$config\" \\
            --nFastqs=$nFastqs \\
            --summary=\"$Summary\" \\
            \"$inputR2\" \"$inputR1\" \\
            -t $SplitcodeNum \\
            -o /dev/null,/dev/null \\
            --keep-grp=\"$SplitcodeLib\" \\
            --keep-r1-r2 --assign \\
            --mapping=\"$Mapping\"
        echo \"splitcode运行完成!\"

        # 5. 等待所有STAR进程完成
        echo \"等待STAR进程完成... (后台进程数: \$(jobs -p | wc -l))\"
        if ! wait; then
            echo \"警告: 部分STAR进程异常退出，检查日志!\"
        fi

        # 6. 卸载基因组
        echo \"卸载基因组...\"
        STAR --genomeDir \"$STARindex\" --genomeLoad Remove --outFileNamePrefix log/GenomeRemove
        
        # 清理FIFO文件
        for s in \"\${samples_arr[@]}\"; do
            rm -f 01-Fastq_Modify_ParaFly/\"lib\${s}_R1.fastq\" 01-Fastq_Modify_ParaFly/\"lib\${s}_R2.fastq\"
        done
        echo \"FIFO文件已清理\"
    "

    echo "===$(date '+%F %T')=== 所有任务完成 ==="
} > "$LogFile" 2>&1