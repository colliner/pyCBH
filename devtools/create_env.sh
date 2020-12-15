#!/usr/bin/env bash

#
# Usage: bash create_env.sh
#

echo 'Creating conda environment for pyCBH'
conda env create -q -f devtools/environment.yml --force
git clone https://github.com/colliner/xyz2mol.git
echo '# To activate this environment, use'
echo '#'
echo '#     $ conda activate pyCBH'
echo '#'
