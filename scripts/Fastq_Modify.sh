rm -rf ./Sample
path=./0-data
files=$(ls $path | grep "\\_R1.fastq.gz")
for filename in $files
do
   echo $filename | awk -F"_R1.fastq.gz" '{print $1;exit}'>> ./Sample
done
wait

mkdir -p ParaFly/
mkdir -p 01-Fastq_Modify_ParaFly
mkdir -p log/
rm -rf ParaFly/01-*.sh
for Sample in $(cat Sample)
do
# Fastq_Modify.sh脚本处理：
# 第一：R1 的CB1有9bp和10bp，依据CB1后紧跟No.10~15bp L6=CAGAGC为判断，符合则CB1补一个C碱基，构成CB1均是10bp;
# 第二：R1和R2的文库序列剪切合并成新的R1，如文库：R1=C10L6U8C10T10  R2=C8L18，则R1剪切No.1~44bp,R2剪切No.1-26；
echo -e \
"
python ${6}/Fastq_Modify.py \
--input_R1 0-data/${Sample}_R1.fastq.gz \
--input_R2 0-data/${Sample}_R2.fastq.gz \
--link_location ${1} \
--link_seq ${2} \
--cut_R1_length ${3} \
--cut_R2_length ${4} \
--output_R1 01-Fastq_Modify_ParaFly/${Sample}_R1.fastq.gz \
--output_R2 01-Fastq_Modify_ParaFly/${Sample}_R2.fastq.gz \
--batch_size  1000 \n">>\
ParaFly/01-Fastq_Modify_ParaFly.sh
done
wait

date
echo 'Start Step1.Fastq_Modfy python!'
for i in 01-Fastq_Modify_ParaFly
do
#-CPU 设定平行任务数，过多平行会减慢单个任务效率
nohup \
ParaFly \
-c ParaFly/${i}.sh \
-CPU ${5} \
-failed_cmds ParaFly/${i}.failed 
done
wait
date
echo 'Step1.Fastq_Modfy python is finished !'

