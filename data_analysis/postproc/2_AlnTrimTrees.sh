#!/bin/bash

###basic calling
#### 2_AlnTrimTrees.sh <gene_list> <seq_root> <species_tree> <threads> <task_id> <out_root>
##### but see the 2_launcher script for more intuiitive startup

###controllers for 2_launcher.sh script/see details there
GENE_LIST=$1
SEQ_ROOT=$2
SPECIES_TREE=$3
THREADS=$4
TASK_ID=$5
OUT_ROOT=${6:-results/gene_centric}

###initi outdirs
mkdir -p "$OUT_ROOT"/{logs,per_gene,alignments_protein,alignments_codon,alignments_trimmed,trees,hyphy,paml}

GENE=$(sed -n "${TASK_ID}p" "$GENE_LIST")
[ -n "$GENE" ] || { echo "No gene for task $TASK_ID"; exit 0; }

GDIR="$OUT_ROOT/per_gene/$GENE"
mkdir -p "$GDIR"
log="$OUT_ROOT/logs/${GENE}.log"
echo "[$(date)] START $GENE" >> "$log"

####collect all prot/cds files, altho protein calls can probably be ditched
#### left those in for posterity although ultimately using CDS to keep 
### QC a little more straightfwd
PROT_FILES=()
CDS_FILES=()

if [ -d "$SEQ_ROOT/$GENE" ]; then
    for f in "$SEQ_ROOT/$GENE"/*; do
        [[ $f == *.fasta || $f == *.fa || $f == *.faa ]] && [[ $f == *protein* ]] && PROT_FILES+=("$f")
        [[ $f == *.fasta || $f == *.fa || $f == *.faa ]] && [[ $f == *cds* || $f == *cdna* ]] && CDS_FILES+=("$f")
    done
else
    for f in "$SEQ_ROOT/${GENE}"*; do
        [[ $f == *.fasta || $f == *.fa ]] && PROT_FILES+=("$f")
    done
fi

[ ${#PROT_FILES[@]} -gt 0 ] || { echo "[$(date)] $GENE: no protein files found, skipping" >> "$log"; exit 1; }

#####merge all cds seqs from primates into one fasta
#
PROT_OUT="$GDIR/original_proteins.fasta"
CDS_OUT="$GDIR/original_cds.fasta"

> "$PROT_OUT"
for f in "${PROT_FILES[@]}"; do
    cat "$f" >> "$PROT_OUT"
done

if [ ${#CDS_FILES[@]} -gt 0 ]; then
    > "$CDS_OUT"
    for f in "${CDS_FILES[@]}"; do
        cat "$f" >> "$CDS_OUT"
    done
fi

echo "[$(date)] Finished collecting all sequences for $GENE" >> "$log"


###### header sanitize for later tools/dedupe only true duplicates
SANITIZED_PROT="$GDIR/original_proteins_safe.fasta"
MAP_PROT="$GDIR/header_map_protein.txt"

awk 'BEGIN{FS="\n"; RS=">"; ORS=""}
NR>1{
    header=$1
    seq=""
    for(i=2;i<=NF;i++){seq=seq $i}
    # Only remove *exact duplicates* across the whole input
    key=header"\n"seq
    if(!seen[key]++){
        print ">"header"\n"seq"\n" >> "'"$SANITIZED_PROT"'"
        print ">"header"\n" >> "'"$MAP_PROT"'"
    }
}' "$PROT_OUT"

if [ ${#CDS_FILES[@]} -gt 0 ]; then
    SANITIZED_CDS="$GDIR/original_cds_safe.fasta"
    MAP_CDS="$GDIR/header_map_cds.txt"
    awk 'BEGIN{FS="\n"; RS=">"; ORS=""}
    NR>1{
        header=$1
        seq=""
        for(i=2;i<=NF;i++){seq=seq $i}
        key=header"\n"seq
        if(!seen[key]++){
            print ">"header"\n"seq"\n" >> "'"$SANITIZED_CDS"'"
            print ">"header"\n" >> "'"$MAP_CDS"'"
        }
    }' "$CDS_OUT"
fi

echo "[$(date)] finished optional header sanitization for $GENE" >> "$log"


#########mafft for trying out protein aln
ALN_PROT="$OUT_ROOT/alignments_protein/${GENE}.aln.faa"

#####3a little diligent pref to use prev written safe version if it exists
#### just an attempt to try to avoid dupes/dirty headers
if [ -f "$GDIR/original_proteins_safe.fasta" ]; then
    PROT_FOR_ALIGN="$GDIR/original_proteins_safe.fasta"
else
    PROT_FOR_ALIGN="$GDIR/original_proteins.fasta"
fi

echo "[$(date)] running mafft for $GENE!" >> "$log"

mafft --thread "$THREADS" \
      --auto \
      "$PROT_FOR_ALIGN" \
      > "$ALN_PROT" 2>> "$log"

if [ $? -ne 0 ] || [ ! -s "$ALN_PROT" ]; then
    echo "[$(date)] uh-oh,mafft failed for $GENE :(" >> "$log"
    exit 1
fi

echo "[$(date)] yay-yay!mafft alignment finished for $GENE" >> "$log"

#######codon-aware aln/macse run !
CODON_ALN="$OUT_ROOT/alignments_codon/${GENE}_codon.aln.fasta"
MACSE_AA="$OUT_ROOT/alignments_protein/${GENE}_macse_aa.fa"

if [ -n "${CDS_OUT:-}" ] && [ -s "$CDS_OUT" ] && command -v macse >/dev/null 2>&1; then
    if [ ! -s "$CODON_ALN" ]; then
        echo "[$(date)] macse codon aln for $GENE" >> "$log"
        macse -prog alignSequences \
            -seq "$CDS_OUT" \
            -out_NT "$CODON_ALN" \
            -out_AA "$MACSE_AA" -fs 20 -stop 20 \
            > "$OUT_ROOT/logs/${GENE}_macse.stdout" \
            2> "$OUT_ROOT/logs/${GENE}_macse.stderr" \
        || { echo "macse failed for $GENE...:(" >> "$log"; }
    fi
else
    echo "[$(date)] no cds for $GENE — skipping macse... :(" >> "$log"
fi

####tree trimming vars to make sure trimal/iqt3 are called to run
TRIMAL_BIN=${TRIMAL_BIN:-trimal}
IQTREE_BIN=${IQTREE_BIN:-iqtree3}

######trimal for prot aln/mafft outs, if needed
PROT_ALN_MAFFT="$OUT_ROOT/alignments_protein/${GENE}.aln.faa"
TRIM_MAFFT="$OUT_ROOT/alignments_trimmed/${GENE}_mafft.trim.faa"

if [ -s "$PROT_ALN_MAFFT" ] && [ ! -s "$TRIM_MAFFT" ]; then
    echo "[$(date)] trimal mafft aln for $GENE" >> "$log"
    "$TRIMAL_BIN" -in "$PROT_ALN_MAFFT" -out "$TRIM_MAFFT" -automated1 > /dev/null 2>> "$OUT_ROOT/logs/${GENE}_trimal_mafft.err" || \
        { echo "trimal mafft failed for $GENE... :(" >> "$log"; exit 1; }
fi


####trimal for prot aln/macse outs, can dump if just using codon
### 
PROT_ALN_MACSE="$OUT_ROOT/alignments_protein/${GENE}_macse_aa.fasta"
TRIM_MACSE="$OUT_ROOT/alignments_trimmed/${GENE}_macse.trim.faa"

if [ -s "$PROT_ALN_MACSE" ] && [ ! -s "$TRIM_MACSE" ]; then
    echo "[$(date)] trimal macse aln for $GENE" >> "$log"
    "$TRIMAL_BIN" -in "$PROT_ALN_MACSE" -out "$TRIM_MACSE" -automated1 > /dev/null 2>> "$OUT_ROOT/logs/${GENE}_trimal_macse.err" || \
        { echo "trimal macse failed for $GENE.. :(" >> "$log"; exit 1; }
fi

####trimal for codon aln/macse codon outs
##33definitely using this so keep!
CODON_ALN="$OUT_ROOT/alignments_codon/${GENE}_codon.aln.fasta"
TRIM_CODON="$OUT_ROOT/alignments_trimmed/${GENE}_codon.trim.fasta"

if [ -s "$CODON_ALN" ] && [ ! -s "$TRIM_CODON" ]; then
    echo "[$(date)] trimal codon aln for $GENE" >> "$log"
    "$TRIMAL_BIN" -in "$CODON_ALN" -out "$TRIM_CODON" -automated1 > /dev/null 2>> "$OUT_ROOT/logs/${GENE}_trimal_codon.err" || \
        { echo "trimal codon failed for $GENE.. :(" >> "$log"; exit 1; }
fi

##iqtree cmd for prot aln trimming
TREE_PREF_MAFFT="$OUT_ROOT/trees/${GENE}_mafft"
if [ -s "$TRIM_MAFFT" ] && [ ! -f "${TREE_PREF_MAFFT}.treefile" ]; then
    echo "[$(date)]iqtree  for mafft-trimmed aln of $GENE" >> "$log"
    "$IQTREE_BIN" -s "$TRIM_MAFFT" -nt AUTO -m MFP -bb 1000 -pre "$TREE_PREF_MAFFT" > "$OUT_ROOT/logs/${GENE}_iqtree_mafft.out" 2> "$OUT_ROOT/logs/${GENE}_iqtree_mafft.err" || \
        { echo "iqtree mafft failed for $GENE... :(" >> "$log"; exit 1; }
fi
TREE_PREF_MACSE="$OUT_ROOT/trees/${GENE}_macse"
if [ -s "$TRIM_MACSE" ] && [ ! -f "${TREE_PREF_MACSE}.treefile" ]; then
    echo "[$(date)] iqtree for macse-trimmed aln of $GENE" >> "$log"
    "$IQTREE_BIN" -s "$TRIM_MACSE" -nt AUTO -m MFP -bb 1000 -pre "$TREE_PREF_MACSE" > "$OUT_ROOT/logs/${GENE}_iqtree_macse.out" 2> "$OUT_ROOT/logs/${GENE}_iqtree_macse.err" || \
        { echo "iqtree macse failed for $GENE... :(" >> "$log"; exit 1; }
fi

###iqtree run for codon aln trimming
######used in draft paper so def keep
TREE_PREF_CODON="$OUT_ROOT/trees/${GENE}_codon"
if [ -s "$TRIM_CODON" ] && [ ! -f "${TREE_PREF_CODON}.treefile" ]; then
    echo "[$(date)] iqtree for codon-trimmed aln of $GENE" >> "$log"
    "$IQTREE_BIN" -s "$TRIM_CODON" -nt AUTO -m MFP -bb 1000 -pre "$TREE_PREF_CODON" > "$OUT_ROOT/logs/${GENE}_iqtree_codon.out" 2> "$OUT_ROOT/logs/${GENE}_iqtree_codon.err" || \
        { echo "iqtree codon failed for $GENE... :(" >> "$log"; exit 1; }
fi
