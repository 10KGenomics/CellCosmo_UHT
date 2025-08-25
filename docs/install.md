# Installation Instructions

## Step 1: Create a conda environment, Install dependent third-party software

```bash
wget https://github.com/10KGenomics/CellCosmo_UHT/blob/main/conda_packages_list.txt
# conda_packages_list.txt file
conda-forge::python=3.10.8
bioconda::star=2.7.10b
bioconda::parafly=r2013_01_21
conda-forge::pigz=2.8
conda-forge::pandas=2.2.3
anaconda::scipy=1.15.3
conda-forge::scanpy=1.11.2
anaconda::openpyxl=3.1.5
conda-forge::jinja2=3.1.6
jmcmurray::os=0.1.4
anaconda::reportlab=3.5.67
conda-forge::svglib=1.5.1
conda-forge::cairosvg=2.8.2     

ENV_NAME=CellCosmo_UHT
conda create -n $ENV_NAME -y --file conda_packages_list.txt

# You can also use mamba installation environment
conda install mamba -y
mamba create -n $ENV_NAME -y --file conda_packages_list.txt
```

## Step 2: Installing CellCosmo_UHT software
```bash
# Download the release package
# Packages under temporary update, Unavailable for the time being
# wget https://github.com/10KGenomics/CellCOSMO_UHT/releases/download/v1.1.1/cell_cosmo_uht-1.0.0.tar.gz
# pip install cell_cosmo_uht-1.0.0.tar.gz

# The current operational approaches are as follows:
# Create a directory named "CellCosmo_UHT" (if it doesn't exist)
mkdir -p CellCosmo_UHT/

# Change the working directory to "CellCosmo_UHT"
cd CellCosmo_UHT/

# Clone the CellCosmo_UHT project from the GitHub repository
git clone https://github.com/10KGenomics/CellCosmo_UHT.git

# Activate the Conda environment (replace $ENV_NAME with the actual environment name)
conda activate $ENV_NAME

# splitcode设置目标安装路径,CellCosmo_UHT环境的bin中
TARGET_DIR="${CellCosmo_UHT_conda_path}/CellCosmo_UHT"

# 下载和解压 splitcode
wget https://github.com/pachterlab/splitcode/archive/refs/tags/v0.31.3.tar.gz
tar -zxvf v0.31.3.tar.gz
cd splitcode-0.31.3/

# 修改源代码：将 jj==0 改为 jj==1
sed -i 's/jj == 0/jj == 1/' src/ProcessReads.cpp
# 创建构建目录
mkdir build
cd build

# 设置编译环境，gcc-env需要详细安装细节
conda create -n gcc-env -y
conda activate gcc-env
conda install -c conda-forge cmake=3.28.0 -y
conda install -c conda-forge gcc_linux-64 -y
conda install -c conda-forge gxx_linux-64 -y
conda install -c conda-forge zlib -y
export CC=$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-gcc
export CXX=$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-g++
export CFLAGS="-D_GNU_SOURCE"
export CXXFLAGS="-D_GNU_SOURCE"

# 配置 CMake，指定安装路径
cmake -DCMAKE_INSTALL_PREFIX="$TARGET_DIR" ..

# 编译和安装
make
make install

# 验证安装
if [ -f "$TARGET_DIR/bin/splitcode" ]; then
    echo "安装成功！splitcode 已安装到: $TARGET_DIR/bin/splitcode"
else
    echo "安装失败，请检查错误信息"
    exit 1
fi

# Run the CellCosmo_UHT script with the `-h` flag to display help documentation
# This shows available parameters and usage instructions for the pipeline
cd ../../
python CellCosmo_UHT/CellCosmo_UHT.py -h
```

