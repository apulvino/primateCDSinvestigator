#!/bin/bash

#SBATCH --account=b1057
#SBATCH --partition=b1057
#SBATCH --job-name=mafft_pal2nal
#SBATCH --output=logs/mafft_pal2nal_%j.log
#SBATCH --error=logs/mafft_pal2nal_%j.err
#SBATCH --time=15:00:00
#SBATCH --mem=32G
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END
#SBATCH --mail-user=anthony.pulvino@northwestern.edu

source activate ortho_explorer

python 2_phyloWork.py
