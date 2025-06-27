# CellCosmo_UHT
Ultra-High-Throughput Single-Cell Transcriptomics Technology
* Propose
   * A collection of bioinfomatics analysis pipelines to process single cell sequencing data generated with Ultra-High-Throughput single-cell RNA products.
* Language
   * Python3(>=3.10.*)
   * R scripts
* Hardware/Software environment
   * x86-64 compatible processors.
   * require at least 30GB of RAM and 4 CPU.
   * centos 7.x 64-bit operating system (Linux kernel 3.10.0, compatible with higher software and hardware configuration).
# Installation
Installation tutorial manual [here](docs/install.md)

# Start Running CellCosmo_UHT
## 1. Build References  genomeDir for Homo sapiens or Mus musculus
Build References genomeDir tutorial manual [here](docs/Build_References_genomeDir.md)

## 2.Run pipeline
```bash
conda activate CellCosmo_UHT
script_path=${path}/CellCosmo_UHT
STARindex=GRCh38_index/
python ${script_path}/CellCosmo_UHT.py \
--script_path ${script_path} \
--SampleName 'nfbmb' \
--SplitLibrary 'False' \
--STARindex ${STARindex} \
--TopCells 8000 \
--Threads 16 \
--CB3_Num 1-20 \
--STARsolo_param '--outReadsUnmapped Fastx --outSAMunmapped Within ' \
--splitCB "1-192" \
--splitSample "nfbmb" \
--splitSample_EstimatedCell_list "3-EstimatedCell_nfbmb.list" \
--splitSample_EstimatedCellMatrix "nfbmb_filtered_feature_bc_matrix"
```

# Support
The officially supported release binaries are available at: (http://www.10kgenomics.com/)
