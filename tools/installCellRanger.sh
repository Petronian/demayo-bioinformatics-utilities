#!/bin/env bash

# Move to home directory.
cd ~

# Download CellRanger.
cellRangerTarGzLink=$(curl -s https://www.10xgenomics.com/support/software/cell-ranger/downloads/previous-versions | sed -nr "s/^.*(https[^\"',]+?cellranger-8\.0\.1\.tar\.gz[^\"',]+).*$/\1/p" | sed -nr 's/\\u0026/\&/gp')
echo "Downloading Cell Ranger 8.0.1..."
curl -o cellranger-8.0.1.tar.gz $cellRangerTarGzLink

# Extract and install CellRanger.
echo "Extracting Cell Ranger 8.0.1..."
tar -xzvf cellranger-8.0.1.tar.gz
rm cellranger-8.0.1.tar.gz

# Add to PATH.
echo "Adding Cell Ranger 8.0.1 to PATH..."
echo "export PATH=~/cellranger-8.0.1/:\$PATH" >> ~/.bashrc

# Download the reference data.
echo "Downloading Cell Ranger human (GRCh38) and mouse (mm10) references..."
curl -o refdata-gex-GRCh38-2020-A.tar.gz https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
curl -o refdata-gex-mm10-2020-A.tar.gz https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz

# Extract the reference data.
echo "Extracting Cell Ranger human (GRCh38) and mouse (mm10) references..."
tar -xvzf refdata-gex-GRCh38-2020-A.tar.gz 
tar -xvzf refdata-gex-mm10-2020-A.tar.gz
rm refdata-gex-GRCh38-2020-A.tar.gz 
rm refdata-gex-mm10-2020-A.tar.gz

echo "Done. All folders located in your home directory (~)."