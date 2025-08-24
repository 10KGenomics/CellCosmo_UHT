# 默认输出文件名
CB3_lib="10K_lib"

# 显示帮助信息
show_help() {
    echo "用法: $0 [选项]"
    echo "根据指定的条形码3编号生成文库文件"
    echo
    echo "选项:"
    echo "  --CB3_Num <编号>    指定编号范围或单独编号（必需）"
    echo "                      示例: '1-8' 或 '1-20 26 28'"
    echo "  --outCB3_lib <文件> 指定输出文件名（默认: 10K_lib）"
    echo "  -h, --help         显示此帮助信息"
    echo
    echo "示例:"
    echo "  $0 --CB3_Num '1-20 26 28' --outCB3_lib my_lib"
}

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case "$1" in
        --CB3_Num)
            CB3_NUM="$2"
            shift 2  # 跳过参数名和值
            ;;
        --outCB3_lib)
            CB3_lib="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "错误: 未知选项 $1"
            show_help
            exit 1
            ;;
    esac
done

# 检查必需参数
if [[ -z "$CB3_NUM" ]]; then
    echo "错误: 必须提供 --CB3_Num 参数"
    show_help
    exit 1
fi

# 清空或创建输出文件
Work_Path=$(pwd)
mkdir -p ${Work_Path}/log
> "${Work_Path}/log/${CB3_lib}"

# 处理编号参数
for arg in $CB3_NUM; do  # 注意这里不加引号以进行单词分割
    if [[ $arg == *-* ]]; then
        # 处理范围格式 (如 1-8)
        IFS='-' read -r start end <<< "$arg"
        
        # 验证是否为数字
        if ! [[ "$start" =~ ^[0-9]+$ && "$end" =~ ^[0-9]+$ ]]; then
            echo "警告: 无效的范围格式 '$arg'，已跳过"
            continue
        fi
        
        # 检查范围有效性
        if (( start > end )); then
            echo "警告: 范围 '$arg' 起始编号大于结束编号，已跳过"
            continue
        fi
        
        # 生成范围内的所有编号
        for ((i = start; i <= end; i++)); do
            echo -e "barcode3_${i},barcode1,linker,barcode2\t01-Fastq_Modify_ParaFly/lib${i}" >> "${Work_Path}/log/${CB3_lib}"
        done
    else
        # 处理单个编号
        if ! [[ "$arg" =~ ^[0-9]+$ ]]; then
            echo "警告: 无效的编号 '$arg'，已跳过"
            continue
        fi
        echo -e "barcode3_${arg},barcode1,linker,barcode2\t01-Fastq_Modify_ParaFly/lib${arg}" >> "${Work_Path}/log/${CB3_lib}"
    fi
done

# 排序处理：按第二列的数字部分排序
# 1. 提取第二列的数字部分
# 2. 按数值排序
# 3. 移除临时添加的排序字段
awk -F'\t' '{ split($2, a, "lib"); print a[2] "\t" $0 }' "${Work_Path}/log/${CB3_lib}" | \
sort -n | \
cut -f2- > "${CB3_lib}.sorted"

# 替换原始文件
mv "${CB3_lib}.sorted" "${Work_Path}/log/${CB3_lib}"

echo "文件已生成：${Work_Path}/log/${CB3_lib}"