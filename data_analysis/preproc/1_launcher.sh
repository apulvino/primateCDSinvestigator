#!/bin/bash

#SBATCH --account=b1057
#SBATCH --partition=b1057
#SBATCH --job-name=getSeqs
#SBATCH --output=getSeqs.out
#SBATCH --error=getSeqs.err
#SBATCH --time=15:00:00
#SBATCH --mem=4G
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END
#SBATCH --mail-user=anthony.pulvino@northwestern.edu

source activate ortho_explorer

#python grabSeqs2.py --config config.yml
python 1_grabSeqs.py
