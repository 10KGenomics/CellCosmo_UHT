# 参考基因组构建STARindex，以GRCh38为例
## STARindex=${Path}/GRCh38_index
```bash
mkdir -p GRCh38_index
cd GRCh38_index

wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar -xzvf refdata-gex-GRCh38-2020-A.tar.gz

STAR --runThreadN 16 \
--runMode genomeGenerate \
--genomeDir ./ \
--genomeFastaFiles refdata-gex-GRCh38-2020-A/fasta/genome.fa \
--sjdbGTFfile      refdata-gex-GRCh38-2020-A/genes/genes.gtf \
--sjdbOverhang 100
```

# 如已安装CellCosmo软件，则运行如下指令
```bash
CellCosmo rna mkref \
--gtf refdata-gex-GRCh38-2020-A/genes/genes.gtf \
--genome-name GRCh38 \
--fasta refdata-gex-GRCh38-2020-A/fasta/genome.fa \
--gene-name-as-name2
```
