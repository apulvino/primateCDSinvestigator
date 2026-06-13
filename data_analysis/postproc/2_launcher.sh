#!/bin/bash

#SBATCH --account=b1057
#SBATCH --partition=b1057
#SBATCH --job-name=gene_cent
#SBATCH --output=results/gene_centric/logs/gene_array_%A_%a.out
#SBATCH --error=results/gene_centric/logs/gene_array_%A_%a.err
#SBATCH --time=13-13:00:00
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END
#SBATCH --array=1-540
#SBATCH --mail-user=anthony.pulvino@northwestern.edu

source activate ortho_explorer

###defining 2_AlnTrimTrees.py vars for launch below/see sister py script
WORKDIR=${WORKDIR:-$(pwd)}
GENELIST=${GENELIST:-gene_list.txt}
SEQ_ROOT=${SEQ_ROOT:-../OrthologExplorer/results_sequences}
SPECIES_TREE=${SPECIES_TREE:-species_tree.nwk}
THREADS=${THREADS:-4}
OUT_ROOT=${OUT_ROOT:-$WORKDIR/results/gene_centric}

###dbl check outdirs init'd, altho should be done in sister pyscript
mkdir -p "$OUT_ROOT"/{logs,per_gene,alignments_protein,alignments_codon,alignments_trimmed,trees,hyphy,paml}

#### a bit of echo to indicate all running smooth
echo "[$(date)] launching gene-centric array jobs..."
echo "gene list: $GENELIST"
echo "sequence root: $SEQ_ROOT"
echo "species tree: $SPECIES_TREE"
echo "threads per task: $THREADS"

####final call to launch align/tree construction/trimming pipeline
bash "$WORKDIR/2_AlnTrimTrees.sh" "$GENELIST" "$SEQ_ROOT" "$SPECIES_TREE" "$THREADS" "$SLURM_ARRAY_TASK_ID" "$OUT_ROOT"
