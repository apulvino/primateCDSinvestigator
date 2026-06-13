#!/usr/bin/env python3

####Bbuilding KEGG-based gene sets in a  clean CSV output, here.
#### ATPulvino 

###load libs
import requests
import time
import csv
from collections import defaultdict

#########3define outputput file and opts for stable extraction of genes
OUTFILE = "gene_sets_union.csv"
SLEEP = 0.35
RETRIES = 3
RETRY_DELAY = 1.5

#######3 define pathways by labeled group and load up the associated kegg ids
PATHWAYS = {
    "bile_acid": {
        "label": "Bile acid biosynthesis & metabolism",
        "kegg": ["hsa00120","hsa00121","hsa04976","hsa04977","hsa04979"]
    },
    "taurine_cysteine": {
        "label": "Taurine and cysteine metabolism",
        "kegg": ["hsa00430","hsa00270","hsa00280", "hsa00920"]
    },
    "glycine": {
        "label": "Glycine metabolism",
        "kegg": ["hsa00260","hsa00670"]
    },
    "retinoid": {
        "label": "Retinoid / vitamin A metabolism",
        "kegg": ["hsa00830","hsa00460"]
    },
    "fatty_cholesterol": {
        "label": "Fatty acid & cholesterol metabolism",
        "kegg": ["hsa01212","hsa04979","hsa00100","hsa00110","hsa01040","hsa00564"]
    }
}

#########defineing a ffew helpers
## fetch url and ensure a couple of retries if connection breaks
def fetch_url(url, retries=RETRIES):
    """Fetch a URL with retries."""
    for attempt in range(retries):
        try:
            r = requests.get(url, timeout=20)
            if r.status_code == 200:
                return r.text
            if r.status_code in (429, 500, 502, 503):
                time.sleep(RETRY_DELAY * (attempt + 1))
        except Exception:
            time.sleep(RETRY_DELAY * (attempt + 1))
    return None

def sleep_polite():
    time.sleep(SLEEP)

#######3pull out the kegg genes for each pathway label/set
def kegg_genes_from_pathway(kegg_path_id):
    """
    Returns a set of (symbol, description) tuples from a KEGG pathway id
    description = name line
    symbol = first element in symbol line
    """
    base = "https://rest.kegg.jp"
    result = set()
    text = fetch_url(f"{base}/link/hsa/{kegg_path_id}")
    if not text:
        return result

    gene_ids = [line.split("\t")[1].split(":")[1] for line in text.strip().splitlines() if "\t" in line]

    for gid in set(gene_ids):
        sleep_polite()
        gene_text = fetch_url(f"{base}/get/hsa:{gid}")
        if not gene_text:
            continue

        symbol, description = None, None
        for line in gene_text.splitlines():
            line = line.strip()
            if line.startswith("SYMBOL") and not symbol:
                ###take 1st symbol
                parts = line.split("SYMBOL")[1].strip()
                symbol = parts.split(",")[0].strip()
            elif line.startswith("NAME") and not description:
                ###name line is description
                parts = line.split("NAME")[1].strip()
                description = parts
            if symbol and description:
                break

        if symbol:
            result.add((symbol, description))
    return result

#####main function for building union
def build_union():
    pathway_to_genes = defaultdict(set)

    for pkey, info in PATHWAYS.items():
        label = info.get("label", pkey)
        print(f"\n== Processing pathway: {pkey} ({label}) ==")
        for k in info.get("kegg", []):
            print(f" KEGG {k} ...", end="", flush=True)
            kgenes = kegg_genes_from_pathway(k)
            sleep_polite()
            print(f" {len(kgenes)} genes")
            pathway_to_genes[pkey].update(kgenes)

    #organzie the csv rows
    rows = []
    for pkey, genes in pathway_to_genes.items():
        label = PATHWAYS[pkey].get("label", pkey)
        for symbol, description in genes:
            rows.append([pkey, label, description, symbol, "KEGG"])

    ###write csv with quotes
    with open(OUTFILE, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pathway_key", "pathway_label", "description", "symbol", "source"])
        for r in rows:
            w.writerow(r)

    print(f"\nWrote {len(rows)} rows to {OUTFILE}")
    return rows

###callmain function/run it
if __name__ == "__main__":
    build_union()
