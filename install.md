# 创建基础环境
conda create -n CellCosmo_UHT python=3.10.8 -y 
# 激活环境
conda activate CellCosmo_UHT 
# 安装生物信息学工具、数据分析库、文件处理工具
conda install -c bioconda star=2.7.10b -y 
conda install bioconda::parafly -y 
conda install conda-forge::pigz -y 
conda install pandas -y 
conda install anaconda::scipy -y 
conda install conda-forge::scanpy -y 
conda install anaconda::openpyxl -y 
