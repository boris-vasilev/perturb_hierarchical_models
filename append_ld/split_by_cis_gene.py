import polars as pl
import os
import math

# ---------------------------
# SETTINGS
# ---------------------------
INPUT = "/rds/project/rds-csoP2nj6Y6Y/biv22/data/pairs/cis_gene_snps.csv"
OUTDIR = "gene_chunks"
GENES_PER_CHUNK = 20
# ---------------------------

df = pl.read_csv(INPUT)

os.makedirs(OUTDIR, exist_ok=True)

# Get unique genes
genes = df.select("cis_gene").unique().to_series().to_list()
genes.sort()

n_chunks = math.ceil(len(genes) / GENES_PER_CHUNK)

print(f"{len(genes)} genes → {n_chunks} chunks")

for i in range(n_chunks):

    start = i * GENES_PER_CHUNK
    end = start + GENES_PER_CHUNK
    gene_subset = genes[start:end]

    chunk_df = df.filter(pl.col("cis_gene").is_in(gene_subset))

    outfile = f"{OUTDIR}/chunk_{i:04d}.csv"
    chunk_df.write_csv(outfile)

    print(f"Wrote {outfile} ({len(gene_subset)} genes)")
