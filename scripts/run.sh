#!/bin/bash

set -e 
set -o pipefail

LOG_FILE="run_log.txt"

{
  echo "===$(date '+%F %T')=== 开始运行 ==="

  /usr/bin/time -v bash -c '
    # 1. Load genome to shared memory
    STAR --genomeDir /home/rs1/2-reference_genome/mm10/STAR_CellCosmo  --genomeLoad LoadAndExit

    # 2. 创建FIFO管道
    for i in {1..21}; do
        mkfifo lib${i}_R1.fastq
        mkfifo lib${i}_R2.fastq
    done

    # 3. 后台并行启动8个 STAR 进程
    for i in {1..21}; do
      STAR --runThreadN 3 \
        --genomeDir /home/rs1/2-reference_genome/mm10/STAR_CellCosmo \
        --outFileNamePrefix star_outs_lib${i}/ \
        --readFilesIn lib${i}_R1.fastq lib${i}_R2.fastq \
        --soloUMIlen 8 \
        --soloCBmatchWLtype 1MM \
        --outSAMattributes CR UR \
        --soloCBwhitelist None \
        --soloFeatures GeneFull GeneFull_Ex50pAS Velocyto \
        --outSAMtype BAM SortedByCoordinate \
        --soloType CB_UMI_Simple \
        --soloCellReadStats Standard \
        --soloCellFilter TopCells 2000 0.99 10 45000 90000 500 0.01 20000 0.001 10000 \
        --outFilterScoreMinOverLread 0.33 \
        --outFilterMatchNminOverLread 0.33 \
        --genomeLoad LoadAndKeep \
        --limitBAMsortRAM 1001013131 &
    done

    # 4. splitcode写入fifo
    splitcode -c /home/jtCai/Scripts/10K-UHT-scRNA-V2.0/config_Demultiplexing.txt \
      --nFastqs=2 \
      --summary=summary.txt \
      ./0-data/JS21_R2.fastq.gz ./0-data/JS21_R1.fastq.gz \
      -t 21 \
      -o /dev/null,/dev/null \
      --keep-grp=10K_lib \
      --keep-r1-r2 --assign \
      --mapping=mapping.txt

    wait

    # 5. 卸载 genome 内存
    STAR --genomeDir /home/rs1/2-reference_genome/mm10/STAR_CellCosmo --genomeLoad Remove
  '

  echo "===$(date '+%F %T')=== 所有任务完成 ==="
} &> "$LOG_FILE"

