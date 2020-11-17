#!/usr/bin/env bash

#
# Usage: bash setup.sh
#

echo 'Creating conda environment for pyCBH'
conda env create -q -f environment.yml --force
git clone https://github.com/colliner/xyz2mol.git

