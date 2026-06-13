#!/usr/bin/env python3

###load libs
import aiohttp
import asyncio
import csv
import os
import async_timeout
import pandas as pd

### define ensembl rest api url,header contetnt, output dir
SERVER = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json"}
OUTDIR = "results_sequences"
os.makedirs(OUTDIR, exist_ok=True)

### reads genesets union from prev 0/script as table sets genes var for symbols
gene_sets_union = pd.read_csv("gene_sets_union.csv")
GENES = ["{}".format(s) for s in gene_sets_union['symbol'].dropna()]

#### pre-defining species list based on official release ensembl genomes avail
SPECIES = [
    "homo_sapiens", "pan_troglodytes", "pan_paniscus", "gorilla_gorilla",
    "pongo_abelii", "nomascus_leucogenys", "chlorocebus_sabaeus", "macaca_mulatta",
    "macaca_fascicularis", "macaca_nemestrina", "papio_anubis", "callithrix_jacchus",
    "saimiri_boliviensis_boliviensis", "aotus_nancymaae", "microcebus_murinus",
    "colobus_angolensis_palliatus", "cercocebus_atys", "theropithecus_gelada",
    "prolemur_simus", "carlito_syrichta", "propithecus_coquereli"
]


manifest_file = os.path.join(OUTDIR, "manifest.csv")
if not os.path.exists(manifest_file):
    with open(manifest_file, "w", newline="") as mf:
        writer = csv.writer(mf)
        writer.writerow(["gene", "species", "transcript_id", "protein_id",
                         "status", "cds_file", "protein_file", "fetch_status"])


#####defining helpers to fetch cds associated meta-attrs for 
##### the genes pulled
async def fetch_json(session, url, retries=3, delay=1.5):
    """retry  json fetch"""
    for i in range(retries):
        try:
            async with async_timeout.timeout(25):
                async with session.get(url, headers=HEADERS) as r:
                    if r.status == 200:
                        return await r.json()
                    elif r.status in [429, 500, 502, 503]:
                        await asyncio.sleep(delay * (i + 1))
        except Exception:
            await asyncio.sleep(delay * (i + 1))
    return None


async def get_gene_id(session, gene, species):
    url = f"{SERVER}/xrefs/symbol/{species}/{gene}?object_type=gene"
    data = await fetch_json(session, url)
    if not data:
        ###fallback try ortholog lookup from human
        ortho_url = f"{SERVER}/homology/symbol/homo_sapiens/{gene}?target_species={species}"
        ortho = await fetch_json(session, ortho_url)
        if ortho and "data" in ortho and ortho["data"]:
            try:
                return ortho["data"][0]["homologies"][0]["id"]
            except Exception:
                return None
        return None
    return data[0]["id"]


async def get_transcripts(session, gene_id, species):
    url = f"{SERVER}/lookup/id/{gene_id}?expand=1"
    data = await fetch_json(session, url)
    if not data or "Transcript" not in data:
        return []
    return [t for t in data["Transcript"] if t.get("biotype") == "protein_coding"]


async def get_sequences(session, tid):
    cds_url = f"{SERVER}/sequence/id/{tid}?type=cds"
    prot_url = f"{SERVER}/sequence/id/{tid}?type=protein"
    cds = await fetch_json(session, cds_url)
    prot = await fetch_json(session, prot_url)
    cds_seq = cds.get("seq") if cds else None
    prot_seq = prot.get("seq") if prot else None
    return cds_seq, prot_seq


def write_manifest(row):
    with open(manifest_file, "a", newline="") as mf:
        csv.writer(mf).writerow(row)


async def process_gene_species(session, gene, species, sem):
    async with sem:
        gene_dir = os.path.join(OUTDIR, gene)
        os.makedirs(gene_dir, exist_ok=True)
        print(f"Processing {gene} in {species}...")

        gid = await get_gene_id(session, gene, species)
        if not gid:
            write_manifest([gene, species, "", "", "", "", "", "gene_not_found"])
            return

        transcripts = await get_transcripts(session, gid, species)
        if not transcripts:
            write_manifest([gene, species, "", "", "", "", "", "no_transcripts"])
            return

        for t in transcripts:
            tid = t["id"]
            pid = t.get("Translation", {}).get("id", "NA")
            status = "canonical" if str(t.get("is_canonical", "0")) == "1" else "alternative"

            cds_seq, prot_seq = await get_sequences(session, tid)
            if not cds_seq or not prot_seq:
                write_manifest([gene, species, tid, pid, status, "", "", "fetch_failed"])
                continue

            header = f">{gene}|{species}|{tid}|{pid}|{status}"
            prot_file = os.path.join(gene_dir, f"{gene}_{species}_protein.fasta")
            cds_file = os.path.join(gene_dir, f"{gene}_{species}_cds.fasta")

            with open(prot_file, "a") as pf:
                pf.write(f"{header}\n{prot_seq}\n")
            with open(cds_file, "a") as cf:
                cf.write(f"{header}\n{cds_seq}\n")

            write_manifest([gene, species, tid, pid, status, cds_file, prot_file, "retrieved"])

###3call to run and organize attr output meta info yoinked from ensembl
async def main():
    sem = asyncio.Semaphore(8)  ###no ]more than 8 concurrent calls to ensembl api
    async with aiohttp.ClientSession() as session:
        tasks = []
        for gene in GENES:
            for sp in SPECIES:
                tasks.append(process_gene_species(session, gene, sp, sem))
        await asyncio.gather(*tasks)


if __name__ == "__main__":
    asyncio.run(main())

