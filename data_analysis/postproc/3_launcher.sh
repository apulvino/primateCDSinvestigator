#!/bin/bash

#SBATCH --account=b1057
#SBATCH --partition=b1057
#SBATCH --job-name=hyphy
#SBATCH --output=/projects/b1057/apulvino/Chapter2/v4_orthexplorer/results/gene_centric/logs/sel_array_%A_%a.out
#SBATCH --error=/projects/b1057/apulvino/Chapter2/v4_orthexplorer/results/gene_centric/logs/sel_array_%A_%a.err
#SBATCH --time=3-13:00:00
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END
#SBATCH --mail-user=anthony.pulvino@northwestern.edu
#SBATCH --array=1-540

source activate ortho_explorer

#####path and vars hardcode setup/init; can be modified as opts 
##### but this isnt a pipeline construction/nextflow proj so leaving hardcoded here...
WORKDIR=${WORKDIR:-$(pwd)}

TRIM_PROTEIN_DIR="/projects/b1057/apulvino/Chapter2/v4_orthexplorer/results/gene_centric/alignments_trimmed"
TRIM_CODON_DIR="/projects/b1057/apulvino/Chapter2/v4_orthexplorer/results/gene_centric/alignments_trimmed"
TREE_DIR="/projects/b1057/apulvino/Chapter2/v4_orthexplorer/results/gene_centric/trees"
HYPHY_OUT_DIR="/projects/b1057/apulvino/Chapter2/v4_orthexplorer/results/gene_centric/hyphy"
MACSE_JAR="macse"

##opt test case for only a single set of cds'/gene-dir
#TEST_GENE=${TEST_GENE:-"AACS"}

##call var for intermediate gene_list prev created and task id access var and
## seq for seq in said intermediate gene-list file
GENE_LIST="gene_list.txt"
TASK_ID=${SLURM_ARRAY_TASK_ID}
SEQ=$(sed -n "${TASK_ID}p" "$GENE_LIST")

if [[ -z "$SEQ" ]]; then
    echo "uh-oh! no gene found for array index $TASK_ID.. :("
    exit 1
fi

echo "array job index: $TASK_ID"
echo "gene: $SEQ"

#####run hyphy prep/analysis script
# -----------------------------
bash 3_HyPhy.sh \
    "$TRIM_PROTEIN_DIR" \
    "$TRIM_CODON_DIR" \
    "$TREE_DIR" \
    "$HYPHY_OUT_DIR" \
    "$SEQ" \
    "$MACSE_JAR"

echo "[$(date)] yay-yay! hyphy complete for $SEQ"
