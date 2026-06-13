#!/bin/bash

#SBATCH --account=b1057
#SBATCH --partition=b1057
#SBATCH --job-name=genel
#SBATCH --output=genel.out
#SBATCH --error=genel.err
#SBATCH --time=15:00:00
#SBATCH --mem=4G
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=END
#SBATCH --mail-user=anthony.pulvino@northwestern.edu

source activate ortho_explorer

#### run script pass the directory/seq_root with all the gene/cds directories holding primate seqs/record out-list txtfile
SEQ_ROOT=${1:-../OrthologExplorer/results_sequences}
OUT_LIST=${2:-gene_list.txt}

####...no redirect before writing
: > "$OUT_LIST"

shopt -s nullglob

for d in "$SEQ_ROOT"/*; do
    if [ -d "$d" ]; then
        files=("$d"/*_protein.fasta "$d"/*_protein.fa "$d"/*.faa "$d"/*.fa "$d"/*.fasta)
        if [ ${#files[@]} -gt 0 ]; then
            echo "$(basename "$d")" >> "$OUT_LIST"
        fi
    fi
done

###fallbackf for single fastas in cds dir/SEQ_ROOT
###this is all wrangling from preproc outdirs so run those first
if [ ! -s "$OUT_LIST" ]; then
    for f in "$SEQ_ROOT"/*.{fa,fasta,faa}; do
        if [ -f "$f" ]; then
            echo "$(basename "${f%.*}")" >> "$OUT_LIST"
        fi
    done
fi

echo "Wrote $(wc -l < "$OUT_LIST") entries to $OUT_LIST"
