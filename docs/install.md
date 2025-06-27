# Installation Instructions

## Step 1: Create a conda environment, Install dependent third-party software

```bash
wget [https://github.com/10KGenomics/CellCosmo_UHT/blob/main/conda_packages_list.txt](https://github.com/10KGenomics/CellCosmo_UHT/blob/main/conda_packages_list.txt)
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
# Packages under temporary update
wget https://github.com/10KGenomics/CellCOSMO_UHT/releases/download/v1.1.1/cell_cosmo_uht-1.0.0.tar.gz
pip install cell_cosmo_uht-1.0.0.tar.gz

# The current operational approaches are as follows:
# Create a directory named "CellCosmo_UHT" (if it doesn't exist)
mkdir -p CellCosmo_UHT/

# Change the working directory to "CellCosmo_UHT"
cd CellCosmo_UHT/

# Clone the CellCosmo_UHT project from the GitHub repository
git clone [https://github.com/10KGenomics/CellCosmo_UHT.git](https://github.com/10KGenomics/CellCosmo_UHT.git)

# Activate the Conda environment (replace $ENV_NAME with the actual environment name)
conda activate $ENV_NAME

# Run the CellCosmo_UHT script with the `-h` flag to display help documentation
# This shows available parameters and usage instructions for the pipeline
python CellCosmo_UHT/CellCosmo_UHT.py -h
```

