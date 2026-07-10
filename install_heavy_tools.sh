#!/bin/bash
echo "=========================================================="
echo "Installing Heavy Bioinformatics Binaries (IQ-TREE & PhiPack)"
echo "=========================================================="

# Ensure basic build tools exist
sudo apt-get update
sudo apt-get install -y build-essential wget unzip

BIN_DIR="/usr/local/bin"

echo "----------------------------------------------------------"
echo "1. Installing IQ-TREE (Maximum Likelihood Phylogenetics)"
echo "----------------------------------------------------------"
cd /tmp
wget -O iqtree.tar.gz https://github.com/iqtree/iqtree2/releases/download/v2.2.2.7/iqtree-2.2.2.7-Linux.tar.gz
tar -xzvf iqtree.tar.gz
sudo cp iqtree-2.2.2.7-Linux/bin/iqtree2 /usr/local/bin/iqtree
echo "IQ-TREE installed successfully."

echo "----------------------------------------------------------"
echo "2. Installing PhiPack (Recombination Index)"
echo "----------------------------------------------------------"
cd /tmp
wget -O phipack.tar.gz https://www.maths.otago.ac.nz/~dbryant/software/PhiPack.tar.gz
mkdir phipack
tar -xzvf phipack.tar.gz -C phipack
cd phipack/src
make
sudo cp Phi /usr/local/bin/Phi
echo "PhiPack installed successfully."

echo "=========================================================="
echo "Installation Complete! Your Web App can now use C++ Tools."
echo "=========================================================="
