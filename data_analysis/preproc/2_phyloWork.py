#!/usr/bin/env python3


#### load libs
import os
import subprocess
import glob
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

#####config out dirs for seqs/align res
OUTDIR = "results_sequences"
ALIGNDIR = "results_alignments"
os.makedirs(ALIGNDIR, exist_ok=True)

###define calls to mafft and pal2nal 
MAFFT_BIN = "mafft"
PAL2NAL_BIN = "pal2nal.pl"

####can be adjusted to match slurm cpu alloc config in companion launcher
CPU_TOTAL = 16
MAFFT_THREADS = 16
MAX_WORKERS = CPU_TOTAL // MAFFT_THREADS


########3few helpers for calling main... 
###### honestly not planning to use much of these
###### ultimately only codon alns gonna be used
def run_mafft(protein_fastas, out_aln):
    """run mafft on concat protein seqs!"""
    tmp_concat = out_aln + ".tmp"
    with open(tmp_concat, "w") as outfile:
        for f in protein_fastas:
            with open(f) as infile:
                outfile.write(infile.read())

    cmd = [MAFFT_BIN, "--thread", str(MAFFT_THREADS), "--auto", tmp_concat]
    print(f"[MAFFT] {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("awww... MAFFT failed:", result.stderr)
        return False

    with open(out_aln, "w") as f:
        f.write(result.stdout)
    os.remove(tmp_concat)
    return True

##### not using pal2nal outs in the final, but leaving it in here anyways/posterity inclusion
def run_pal2nal(prot_aln, cds_fastas, out_file):
    """run pal2nal for codon aln back translation"""
    tmp_cds = out_file + ".cds.tmp"
    with open(tmp_cds, "w") as outfile:
        for f in cds_fastas:
            with open(f) as infile:
                outfile.write(infile.read())

    cmd = [PAL2NAL_BIN, prot_aln, tmp_cds, "-output", "fasta"]
    print(f"[PAL2NAL] {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("aww... PAL2NAL failed:", result.stderr)
        os.remove(tmp_cds)
        return False

    with open(out_file, "w") as f:
        f.write(result.stdout)
    os.remove(tmp_cds)
    return True
def process_gene(gene):
    """align one gene across species, return summary message"""
    gene_dir = os.path.join(OUTDIR, gene)
    prot_fastas = sorted(glob.glob(os.path.join(gene_dir, f"{gene}_*_protein.fasta")))
    cds_fastas = sorted(glob.glob(os.path.join(gene_dir, f"{gene}_*_cds.fasta")))

    if len(prot_fastas) < 3:
        return f"Skipping {gene}: not enough sequences ({len(prot_fastas)})"

    out_prot_aln = os.path.join(ALIGNDIR, f"{gene}_protein_aln.fasta")
    out_codon_aln = os.path.join(ALIGNDIR, f"{gene}_codon_aln.fasta")

    print(f"\n=== Processing {gene} ===")
    if run_mafft(prot_fastas, out_prot_aln):
        if run_pal2nal(out_prot_aln, cds_fastas, out_codon_aln):
            return f"yay, {gene} is  done! ({len(prot_fastas)} seqs)"
        else:
            return f"aww... PAL2NAL failed for {gene}"
    else:
        return f"aww... MAFFT failed for {gene}"


#####call to main for running 
if __name__ == "__main__":
    genes = [g for g in os.listdir(OUTDIR) if os.path.isdir(os.path.join(OUTDIR, g))]
    print(f"Detected {len(genes)} gene directories")

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(process_gene, g): g for g in genes}
        for fut in tqdm(as_completed(futures), total=len(futures), desc="Alignments"):
            print(fut.result())
