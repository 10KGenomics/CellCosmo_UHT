
conda create -n CellCosmo_UHT python=3.10.8 
conda activate CellCosmo_UHT
#
conda install -c bioconda star=2.7.10b -y
conda install bioconda::parafly -y
conda install conda-forge::pigz -y 
conda install pandas -y
conda install anaconda::scipy -y 
conda install conda-forge::scanpy -y
conda install anaconda::openpyxl -y