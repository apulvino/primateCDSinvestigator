#!/usr/bin/env python3


###load libs
import gzip
import re
from pathlib import Path
import pandas as pd

####set tree path/output metrics file for data pulled
TREES_DIR = Path("/projects/b1057/apulvino/Chapter2/v4_orthexplorer/results/gene_centric/trees")
OUTPUT_FILE = Path("iqtree_metrics.csv")

def open_file(path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path, "r")

def parse_files(iqtree_path, log_path=None):
    gene_id = iqtree_path.stem.replace("_codon", "")

    row = {
        "gene_id": gene_id,

        #####pullout likelihood info
        "raw_lnL": None,
        "raw_consensus_lnL": None,
        "raw_unconstrained_lnL": None,

        ####pull compledxity
        "raw_n_params": None,

        ####pull out tree metric/discordance info
        "raw_tree_length": None,
        "raw_internal_length": None,
        "rf_distance_ml_vs_consensus": None,
        "n_splits": None,

        ###pullout the ICs
        "raw_AIC": None,
        "raw_AICc": None,
        "raw_BIC": None,

        ####pull aln stats
        "n_sequences": None,
        "n_sites": None,
        "n_constant_sites": None,
        "pct_constant_sites": None,
        "n_invariant_sites": None,
        "n_parsimony_informative_sites": None,
        "n_singleton_sites": None,
        "n_distinct_site_patterns": None,
        "prop_invariable_sites": None,

        ###Pull analysnis note bootstrap info
        "analysis_type": None,
        "ufboot_replicates": None,
        "n_bootstrap_trees": None,

        ####pull QC info
        "n_high_gap_sequences": None,
        "has_high_gap_warning": False,
    }
###painstaking pattern matcher lib
    patterns = {
        "raw_lnL": r"Log-likelihood of the tree:\s+([-\d\.]+)",
        "raw_consensus_lnL": r"Log-likelihood of consensus tree:\s+([-\d\.]+)",
        "raw_unconstrained_lnL": r"Unconstrained log-likelihood.*:\s+([-\d\.]+)",
        "raw_n_params": r"Number of free parameters.*:\s+(\d+)",
        "raw_tree_length": r"Total tree length.*:\s+([-\d\.]+)",
        "raw_internal_length": r"Sum of internal branch lengths.*:\s+([-\d\.]+)",
        "rf_distance_ml_vs_consensus": r"Robinson-Foulds distance.*:\s+(\d+)",
        "raw_AIC": r"Akaike information criterion \(AIC\) score:\s+([-\d\.]+)",
        "raw_AICc": r"Corrected Akaike information criterion \(AICc\) score:\s+([-\d\.]+)",
        "raw_BIC": r"Bayesian information criterion \(BIC\) score:\s+([-\d\.]+)",

        "analysis_type": r"Type of analysis:\s+(.+)",
        "ufboot_replicates": r"ultrafast bootstrap\s+\((\d+)\s+replicates\)",
        "n_bootstrap_trees": r"Consensus tree is constructed from\s+(\d+)\s+bootstrap",

        "n_sequences": r"Input data:\s+(\d+)\s+sequences",
        "n_sites": r"Input data:.*with\s+(\d+)\s+nucleotide sites",
        "n_constant_sites": r"Number of constant sites:\s+(\d+)",
        "pct_constant_sites": r"Number of constant sites:.*=\s+([\d\.]+)%",
        "n_invariant_sites": r"Number of invariant.*sites:\s+(\d+)",
        "n_parsimony_informative_sites": r"Number of parsimony informative sites:\s+(\d+)",
        "n_distinct_site_patterns": r"Number of distinct site patterns:\s+(\d+)",
        "n_splits": r"\d+\s+taxa and\s+(\d+)\s+splits",
    }

    ####walk thru iqtree file forest
    with open_file(iqtree_path) as f:
        for line in f:
            line = line.strip()
            for key, pat in patterns.items():
                if row[key] is None:
                    m = re.search(pat, line)
                    if m:
                        row[key] = m.group(1)

    #####incorporation of informative/invar/gappy-site logging 
    if log_path and log_path.exists():
        with open_file(log_path) as f:
            for line in f:
                line = line.strip()

                if row["n_singleton_sites"] is None:
                    m = re.search(r"\d+\s+parsimony-informative,\s+(\d+)\s+singleton sites", line)
                    if m:
                        row["n_singleton_sites"] = m.group(1)

                if row["prop_invariable_sites"] is None:
                    m = re.search(r"Proportion of invariable sites:\s+([\d\.]+)", line)
                    if m:
                        row["prop_invariable_sites"] = m.group(1)

                if "WARNING:" in line and "gaps/ambiguity" in line:
                    row["has_high_gap_warning"] = True
                    m = re.search(r"WARNING:\s+(\d+)\s+sequences contain", line)
                    if m:
                        row["n_high_gap_sequences"] = m.group(1)

    return row

###### main call,iter over treefiles and send thru main parser fxn
rows = []

for iqtree in sorted(TREES_DIR.glob("*_codon.iqtree*")):
    log = iqtree.with_suffix(".log")
    row = parse_files(iqtree, log if log.exists() else None)
    rows.append(row)

df = pd.DataFrame(rows)
df.to_csv(OUTPUT_FILE, index=False)
print(f"stats written to {OUTPUT_FILE} ! :) ")
