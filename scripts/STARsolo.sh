rm -rf ./Sample
path=./0-data
files=$(ls $path | grep "\\_R1.fastq.gz")
for filename in $files
do
   echo $filename | awk -F"_" '{print $1;exit}'>> ./Sample
done
wait

# 确定物种
STARindex=${1}
# 每个CB3子会计数一个CB3子文库=？细胞数
TopCells=${2}
# STARsolo的线程数
Threads=${3}
# CB3范围区间
CB3_Num=${4}
# STARsolo 平行任务数，默认1，需结合线程数设定
ParallelCPU_Step1=${5}
STARsolo_param=${7}

# 选取CB3的CB list
rm -rf ./Scripts
mkdir -p ./Scripts
cp -r ${6}/Barcode ./Scripts/
rm -rf Scripts/Barcode/CB3
mkdir  Scripts/Barcode/CB3
IFS='-' read start end <<< "$CB3_Num"
lines=$(awk -v s="$start" -v e="$end" 'NR>=s && NR<=e' Scripts/Barcode/CB3.ref.list)
echo "$lines" > Scripts/Barcode/CB3.target.list

for CB in $lines
do
echo -e ${CB} > Scripts/Barcode/CB3/${CB}.list
done
wait


mkdir -p 02-STARsolo
rm -rf ParaFly/02-*.sh
#
if [ -n "${STARsolo_param}" ]; then
# 参数存在，使用STARsolo_param
for Sample in $(cat Sample)
do
for CBID in $(cat Scripts/Barcode/CB3.target.list)
do
echo -e \
"STAR \
--genomeDir ${STARindex} \
--readFilesIn \
01-Fastq_Modify_ParaFly/${Sample}_R2.fastq.gz \
01-Fastq_Modify_ParaFly/${Sample}_R1.fastq.gz \
--soloCBwhitelist \
Scripts/Barcode/CB1.list  \
Scripts/Barcode/CB2.list  \
Scripts/Barcode/CB3/${CBID}.list \
--soloCellFilter  TopCells ${TopCells} \
0.99 10 45000 90000 500 0.01 20000 0.001 10000 \
--soloType CB_UMI_Complex \
--soloCBposition 0_0_0_9 0_24_0_33 0_44_0_51 \
--soloUMIposition 0_16_0_23 \
--soloUMIlen 8  \
--soloCBmatchWLtype 1MM \
--soloCellReadStats Standard \
--soloBarcodeReadLength 0 \
--outFileNamePrefix 02-STARsolo/${Sample}_${CBID}- \
--runThreadN ${Threads} \
--outSAMtype SAM \
--outSAMunmapped Within \
--soloFeatures Gene GeneFull_Ex50pAS Velocyto \
--readFilesCommand zcat \
--outFilterMatchNmin 0 \
--outFilterScoreMinOverLread 0.33 \
--outFilterMatchNminOverLread 0.33 \
${STARsolo_param} \n">>\
ParaFly/02-STARsolo.sh
done
wait
done
wait
else
for Sample in $(cat Sample)
do
for CBID in $(cat Scripts/Barcode/CB3.target.list)
do
echo -e \
"STAR \
--genomeDir ${STARindex} \
--readFilesIn \
01-Fastq_Modify_ParaFly/${Sample}_R2.fastq.gz \
01-Fastq_Modify_ParaFly/${Sample}_R1.fastq.gz \
--soloCBwhitelist \
Scripts/Barcode/CB1.list  \
Scripts/Barcode/CB2.list  \
Scripts/Barcode/CB3/${CBID}.list \
--soloCellFilter  TopCells ${TopCells} \
0.99 10 45000 90000 500 0.01 20000 0.001 10000 \
--soloType CB_UMI_Complex \
--soloCBposition 0_0_0_9 0_24_0_33 0_44_0_51 \
--soloUMIposition 0_16_0_23 \
--soloUMIlen 8  \
--soloCBmatchWLtype 1MM \
--soloCellReadStats Standard \
--soloBarcodeReadLength 0 \
--outFileNamePrefix 02-STARsolo/${Sample}_${CBID}- \
--runThreadN ${Threads} \
--outSAMtype None \
--outSAMunmapped Within \
--soloFeatures Gene GeneFull_Ex50pAS Velocyto \
--readFilesCommand zcat \
--outFilterMatchNmin 0 \
--outFilterScoreMinOverLread 0.33 \
--outFilterMatchNminOverLread 0.33 \n">>\
ParaFly/02-STARsolo.sh
done
wait
done
wait
fi
date

echo 'Start Step2.STARsolo !'
for i in 02-STARsolo
do
#-CPU 设定平行任务数，过多平行会减慢单个任务效率
nohup \
ParaFly \
-c ParaFly/${i}.sh \
-CPU ${ParallelCPU_Step1} \
-failed_cmds ParaFly/${i}.failed 
done
wait
date
echo 'Step2.STARsolo is finished !'