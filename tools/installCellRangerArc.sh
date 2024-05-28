#!/bin/env bash

# Move to home directory.
cd ~

# Download CellRanger.
cellRangerArcTarGzLink=$(curl -s https://www.10xgenomics.com/support/software/cell-ranger-arc/downloads/previous-versions | sed -nr "s/^.*(https[^\"',]+?cellranger-arc-2\.0\.2\.tar\.gz[^\"',]+).*$/\1/p" | sed -nr 's/\\u0026/\&/gp')
echo "Downloading Cell Ranger ARC 2.0.2..."
curl -o cellranger-arc-2.0.2.tar.gz $cellRangerArcTarGzLink

# Extract and install CellRanger.
echo "Extracting Cell Ranger ARC 2.0.2..."
tar -xzvf "cellranger-arc-2.0.2.tar.gz"
rm "cellranger-arc-2.0.2.tar.gz"

# Add to PATH.
echo "Adding Cell Ranger ARC 2.0.2 to PATH..."
echo "export PATH=~/cellranger-arc-2.0.2/:\$PATH" >> ~/.bashrc

# Download the reference data.
echo "Downloading Cell Ranger ARC human (GRCh38) and mouse (GRCm38/mm10) reference..."
curl -o refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz "https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz"
curl -o refdata-cellranger-arc-mm10-2020-A-2.0.0.tar.gz "https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-mm10-2020-A-2.0.0.tar.gz"

# Extract the reference data.
echo "Extracting Cell Ranger ARC human (GRCh38) and mouse (GRCm38/mm10) reference..."
tar -xvzf refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz 
tar -xvzf refdata-cellranger-arc-mm10-2020-A-2.0.0.tar.gz 
rm refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz 
rm refdata-cellranger-arc-mm10-2020-A-2.0.0.tar.gz

echo "Done. All folders located in your home directory (~)."