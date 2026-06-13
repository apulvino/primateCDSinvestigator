#!/bin/bash

####input codon aln info/tree/seqs for hyphy
TRIM_PROTEIN_DIR="$1"      ####like alignments_trimmed/protein
TRIM_CODON_DIR="$2"        #####lijke alignments_trimmed/codon
TREE_DIR="$3"              ####the iqtree out .treefiles
HYPHY_OUT_DIR="$4"         #####hyphy outs location
SEQ="$5"  
MACSE_JAR=${6:-"macse"}

mkdir -p "$HYPHY_OUT_DIR"
HYPHY_ALIGN_DIR="$HYPHY_OUT_DIR/alignments"
mkdir -p "$HYPHY_ALIGN_DIR"

## init tmp dir although not really using this anymore actually
TMP_DIR="$HYPHY_OUT_DIR/tmp"
mkdir -p "$TMP_DIR"

####infile alignments trimmed, trees, fastas for hyphy cmd use and
###hardcoded cpu usage
PROTEIN_FILE="$TRIM_PROTEIN_DIR/${SEQ}_mafft.trim.faa"
CODON_FILE="$TRIM_CODON_DIR/${SEQ}_codon.trim.fasta"
TREE_FILE="$TREE_DIR/${SEQ}_codon.treefile"
HYPHY_CODON="$HYPHY_ALIGN_DIR/${SEQ}_for_hyphy.fasta"
CPU=4

###check input to make sure we're not 'sending in' trash
if [ ! -s "$PROTEIN_FILE" ]; then
    echo "uh-oh! protein aln not found: $PROTEIN_FILE"
    exit 1
fi

if [ ! -s "$CODON_FILE" ]; then
    echo "uh-oh! codon aln not found: $CODON_FILE"
    exit 1
fi

if [ ! -s "$TREE_FILE" ]; then
    echo "uh-oh! treefile not found: $TREE_FILE"
    exit 1
fi

####mask internal stops and frameshifts
###had lots of trouble doing for prot so trudging ahead w/
####just the codon data...
echo "[$(date)] Masking internal stops with MACSE..."
"$MACSE_JAR" -prog exportAlignment \
-align "$CODON_FILE" \
-out_NT "$HYPHY_CODON" \
-codonForInternalStop NNN \
-codonForInternalFS NNN \
-codonForFinalStop NNN \
-gc_def 1

####quick .fa cleanup for picky hyphy tools
sed -i 's/|/_/g' "$HYPHY_CODON"

#####run gard for infer recomb breakpoint
GARD_OUT="$HYPHY_OUT_DIR/${SEQ}_GARD.json"
hyphy gard \
--alignment "$HYPHY_CODON" ENV="TOLERATE_NUMERICAL_ERRORS=1;" \
--tree "$TREE_FILE" \
--output "$GARD_OUT" \
--CPU "$CPU"

GARDIN="${HYPHY_ALIGN_DIR}/${SEQ}_for_hyphy.fasta.best-gard"
######run of fel site level inference
###### assign codon aln and tree paths/use best inferred gard bp model
FEL_OUT="$HYPHY_OUT_DIR/${SEQ}_FEL.json"
echo "running hyphy-fel"
hyphy fel --alignment "$GARDIN" \
--output "$FEL_OUT" \
--CPU "$CPU"

######run of meme site level inference 
MEME_OUT="$HYPHY_OUT_DIR/${SEQ}_MEME.json"
echo "running hyphy-meme"
hyphy meme --alignment "$GARDIN" \
--output "$MEME_OUT" \
--CPU "$CPU"

######run of busted tree-wide inference
BUSTED_OUT="$HYPHY_OUT_DIR/${SEQ}_BUSTED.json"
hyphy busted \
--alignment "$GARDIN" \
--output "$BUSTED_OUT" \
--CPU "$CPU"

######run of fubar tree-wide inference
FUBAR_OUT="$HYPHY_OUT_DIR/${SEQ}_FUBAR.json"
hyphy fubar \
--alignment "$GARDIN" \
--output "$FUBAR_OUT" \
--CPU "$CPU"

###############################################################################################
###############################################################################################
################### aBSREL AND RELAX CANNOT BE RUN ON MULTI-PART GARD INPUT###################
###############################################################################################
################################################################################################
######sadly means you gotta split the gard output for all alignments per partition into output nexus file
python split_gard_nexus_by_partition.py "$GARDIN" \
--gard-json "$GARD_OUT" --outdir "$HYPHY_OUT_DIR"/alignments/parts

nex_files="${HYPHY_OUT_DIR}/alignments/parts/${SEQ}.fasta_partition_"*.nex

for f in $nex_files; do \
	[ -f "$f" ] || continue    ##skip empty glob condish
	part_ext=${f##*partition_}
	part_num=${part_ext%.nex}

	ALIGN="$f"
	aBSREL_OUT="${HYPHY_OUT_DIR}/${SEQ}part${part_num}_aBSREL.json"

	echo "Running absrel for ${SEQ} partition ${part_num} -> ${aBSREL_OUT}"
	hyphy absrel \
	--alignment "$ALIGN" \
	--output "$aBSREL_OUT" \
	--CPU "$CPU" || {
		echo "uh-oh! hyphy failed for ${ALIGN}...:(" >&2
		continue
	}
done

################## WORKING TEMPLATE... definitely not finished...
# ##running relax
## variables you already have per-seq and per-partition
##PART_NEX="parts/${SEQ}_partition_${P}.nex"
#LABELED_NWK="tmp/${SEQ}_part${P}_labeled.nwk"   ###from hyphy label-tree
#RELAX_NEX="relax_inputs/${SEQ}_part${P}_RELAX.nex"
#
#python3 make_relax_nexus_from_partition.py \
#--partition-nex "$PART_NEX" \
#--labeled-tree "$LABELED_NWK" \
#--out "$RELAX_NEX"
#
#hyphy relax \
#--alignment "$RELAX_NEX" \
#--test TEST \
#--reference-group REFERENCE \
#--output "${HYPHY_OUT_DIR}/${SEQ}part${part_num}_RELAX.json" \
#--CPU "$CPU"
