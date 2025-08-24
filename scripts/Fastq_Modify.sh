#!/bin/bash

# 显示帮助信息
show_help() {
    echo "用法: $0 [选项] <CB3列表文件> <CB3范围> <最大错配> <CB3位置> <Link位置> <Link序列> <R1保留长度> <R2保留长度> <并行任务数> <脚本路径>"
    echo
    echo "功能:"
    echo "  本脚本用于并行处理FASTQ文件，根据CB3条形码拆分文件并进行序列修改"
    echo
    echo "参数说明:"
    echo "  <CB3列表文件>    包含CB3条形码的文件路径"
    echo "  <CB3范围>        使用的条形码范围（格式：起始-结束，如1-18）"
    echo "  <最大错配>       允许的最大错配数（0-4）"
    echo "  <CB3位置>        CB3在R2中的位置（格式：起始-结束，如1-8）"
    echo "  <Link位置>       Link序列在R1中的位置（格式：起始-结束，如10-15）"
    echo "  <Link序列>       需要匹配的Link序列（如CAGAGC）"
    echo "  <R1保留长度>     保留R1开头的碱基数"
    echo "  <R2保留长度>     保留R2开头的碱基数"
    echo "  <并行任务数>     并行处理的任务数量"
    echo "  <脚本路径>       Python脚本Fastq_SplitCB3_Modify.py的路径"
    echo
    echo "选项:"
    echo "  -h, --help       显示此帮助信息"
    echo
    echo "示例:"
    echo "  $0 CB3.ref.list 1-18 2 1-8 10-15 CAGAGC 44 26 8 /path/to/scripts"
    echo
    echo "输出:"
    echo "  在01-Fastq_Modify_ParaFly目录中为每个样本创建处理后的FASTQ文件"
    echo
    echo "依赖:"
    echo "  需要安装ParaFly并行处理工具"
    echo
    exit 0
}

# 检查是否请求帮助
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    show_help
fi

# 检查参数数量
if [ $# -ne 10 ]; then
    echo "错误: 需要10个参数，但提供了 $# 个参数"
    echo "使用 -h 查看帮助信息"
    exit 1
fi

# 记录开始时间
start_time=$(date +%s)
echo "脚本开始时间: $(date -d @$start_time '+%Y-%m-%d %H:%M:%S')"

# 设置工作路径
Work_Path=$(pwd)
echo "当前工作目录: $Work_Path"

# 获取参数
CB3list="$1"
CB3_Num="$2"
Mismatch="$3"
CB3_location="$4"
link1_location="$5"
link1_seq="$6"
cut_R1_length="$7"
cut_R2_length="$8"
ParallelCPU_Fastq_Modify="$9"
script_path="${10}"

# 验证参数有效性
echo "验证参数..."
if ! [[ "$Mismatch" =~ ^[0-4]$ ]]; then
    echo "错误: 最大错配数必须为0-4之间的整数"
    exit 1
fi

if ! [[ "$ParallelCPU_Fastq_Modify" =~ ^[0-9]+$ ]]; then
    echo "错误: 并行任务数必须是正整数"
    exit 1
fi

if [ ! -f "$CB3list" ]; then
    echo "错误: CB3列表文件不存在: $CB3list"
    exit 1
fi

if [ ! -d "$script_path" ]; then
    echo "错误: 脚本路径不存在: $script_path"
    exit 1
fi

# 创建必要的目录
mkdir -p {ParaFly,01-Fastq_Modify_ParaFly,log} 2>/dev/null
echo "创建目录: ParaFly, 01-Fastq_Modify_ParaFly, log"

# 清理旧文件
rm -f "${Work_Path}/Sample" ParaFly/01-*.sh 2>/dev/null

# 查找样本
echo "在0-data目录中查找样本..."
path="${Work_Path}/0-data"
if [ ! -d "$path" ]; then
    echo "错误: 数据目录不存在: $path"
    exit 1
fi

# 更安全地处理带空格的文件名
find "$path" -maxdepth 1 -name "*_R1.fastq.gz" -exec basename {} \; | \
    sed 's/_R1\.fastq\.gz//' | \
    sort -u > "${Work_Path}/Sample"

sample_count=$(wc -l < "${Work_Path}/Sample")
if [ "$sample_count" -eq 0 ]; then
    echo "错误: 在$path目录中未找到样本文件"
    echo "请确保文件名格式为: <样本名>_R1.fastq.gz"
    exit 1
fi

echo "找到 $sample_count 个样本"

# 为每个样本创建处理命令
echo "生成并行任务脚本..."
while IFS= read -r Sample; do
    cat <<- EOF >> ParaFly/01-Fastq_Modify_ParaFly.sh
    nohup python "${script_path}/Fastq_SplitCB3_Modify.py" \
	 --Sample "${Sample}" \
	 --CB3list "${CB3list}" \
	 --CB3_Num "${CB3_Num}" \
	 --Mismatch "${Mismatch}" \
	 --Location "${CB3_location}" \
	 --link_location "${link1_location}" \
	 --link_seq "${link1_seq}" \
	 --cut_R1_length "${cut_R1_length}" \
	 --cut_R2_length "${cut_R2_length}" \
	 --batch_size 10000 \
	 --inputR1  "${path}/${Sample}_R1.fastq.gz" \
	 --inputR2  "${path}/${Sample}_R2.fastq.gz" \
	 --outDir   "${Work_Path}/01-Fastq_Modify_ParaFly/${Sample}" \
	 >> "log/${Sample}_Fastq_Modify.log" 2>&1
	EOF
done < "${Work_Path}/Sample"

# 检查是否生成了任务脚本
if [ ! -s "ParaFly/01-Fastq_Modify_ParaFly.sh" ]; then
    echo "错误: 未能生成任务脚本"
    exit 1
fi

echo "开始处理FASTQ文件..."
echo "并行任务数: $ParallelCPU_Fastq_Modify"
echo "详细日志请查看log目录"

# 使用ParaFly并行处理
ParaFly \
    -c ParaFly/01-Fastq_Modify_ParaFly.sh \
    -CPU "$ParallelCPU_Fastq_Modify" \
    -failed_cmds ParaFly/01-Fastq_Modify_ParaFly.failed_cmds \
    -v

# 检查失败命令
if [ -s "ParaFly/01-Fastq_Modify_ParaFly.failed_cmds" ]; then
    failed_count=$(wc -l < "ParaFly/01-Fastq_Modify_ParaFly.failed_cmds")
    echo "警告: $failed_count 个任务失败，请查看 ParaFly/01-Fastq_Modify_ParaFly.failed_cmds"
else
    echo "所有任务成功完成"
fi

# 计算运行时间
end_time=$(date +%s)
duration=$((end_time - start_time))
hours=$((duration / 3600))
minutes=$(( (duration % 3600) / 60 ))
seconds=$((duration % 60))

echo "处理完成!"
echo "总运行时间: ${hours}小时 ${minutes}分钟 ${seconds}秒"
echo "结束时间: $(date -d @$end_time '+%Y-%m-%d %H:%M:%S')"