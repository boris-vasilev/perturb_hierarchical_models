#!/usr/bin/env python
import scanpy as sc
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
from tqdm import tqdm
import mygene


def map_ensg_to_symbol(df, gene_col='effect'):
    """
    Maps ENSEMBL IDs in `gene_col` to gene symbols using mygene.

    Returns a DataFrame with:
        - 'effect_ensg' = original ENSEMBL
        - 'effect' = gene symbol
    """
    mg = mygene.MyGeneInfo()
    ensg_list = df[gene_col].unique().tolist()

    print(f"Mapping {len(ensg_list)} ENSG IDs to symbols...")
    results = mg.querymany(ensg_list, scopes='ensembl.gene', fields='symbol', species='human')

    # Build mapping dictionary
    symbol_map = {r['query']: r.get('symbol', None) for r in results if 'symbol' in r}

    # Rename original column and add mapped symbol
    df = df.copy()
    df = df.rename(columns={gene_col: 'effect_ensg'})
    df['effect'] = df['effect_ensg'].map(symbol_map)

    # Remove rows where symbol mapping failed
    df = df.dropna(subset=['effect'])
    return df

def extract_DEGs(DEG_dir, significant=True, bulk=True):
    DEG_dir = Path(DEG_dir)
    DEG_files = list(DEG_dir.glob("*.tsv"))
    if len(DEG_files) == 0:
        raise FileNotFoundError(f"No DEG TSV files found in {DEG_dir}")
    
    dfs = []
    for f in DEG_files:
        DEGs = pd.read_csv(f, sep="\t", index_col=0)
        perturbed_gene = f.stem
        if bulk:
            # mimic column renaming from R
            DEGs = DEGs.rename(columns={"log2FoldChange": "avg_log2FC",
                                        "padj": "p_val_adj"})
        DEGs = DEGs.reset_index().rename(columns={"gene": "effect"})
        DEGs['perturbation'] = perturbed_gene
        dfs.append(DEGs)
    
    all_DEGs = pd.concat(dfs, ignore_index=True)
    
    if significant:
        all_DEGs = all_DEGs[all_DEGs['p_val_adj'] < 0.05]
    
    # Keep only columns needed for QC
    pairs = all_DEGs[['perturbation', 'effect']]
    return pairs

# -----------------------------
# ARGUMENTS
# -----------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--cells", type=str, required=True, help="Cell type/line")
args = parser.parse_args()
cells = args.cells

# -----------------------------
# FILE PATHS
# -----------------------------
DATA_DIR = Path("/rds/project/rds-csoP2nj6Y6Y/biv22/perturb_pseudobulk_de/data")
PERTURB_DIR = Path(f"/rds/project/rds-csoP2nj6Y6Y/biv22/data/perturb/{cells}")
OUT_DIR = Path(f"/rds/project/rds-csoP2nj6Y6Y/biv22/data/perturb_QC/{cells}")
OUT_DIR.mkdir(exist_ok=True, parents=True)

h5ad_files = {
    "Jurkat": "jurkat_raw_singlecell_01.h5ad",
    "K562_GenomeWide": "K562_gwps_raw_singlecell.h5ad",
    "K562_essential": "K562_essential_raw_singlecell_01.h5ad",
    "HepG2": "hepg2_raw_singlecell_01.h5ad",
    "RPE1": "rpe1_raw_singlecell_01.h5ad"
}

# -----------------------------
# LOAD DATA
# -----------------------------
adata = sc.read_h5ad(DATA_DIR / h5ad_files[cells])

# Ensure perturbation is in obs
if 'gene' not in adata.obs:
    raise ValueError("Expected 'gene' column in adata.obs for perturbations")

# -----------------------------
# LOAD DEGs / extract pairs
# -----------------------------
# Assumes you have an equivalent extract_DEGs output as CSVs in PERTURB_DIR
DEG_dir = f"/rds/project/rds-csoP2nj6Y6Y/biv22/data/perturb/{cells}"
pairs = extract_DEGs(DEG_dir, significant=False, bulk=True)
pairs = map_ensg_to_symbol(pairs)


# -----------------------------
# GENE QC FUNCTION
# -----------------------------
def gene_qc(adata, pairs, expression_threshold=0, min_cells=7):
    # Filter pairs to valid genes and perturbations
    valid_genes = adata.var_names
    valid_perts = adata.obs['gene'].unique()
    pairs = pairs[pairs['effect_ensg'].isin(valid_genes) & pairs['perturbation'].isin(valid_perts)]
    print(f"Total valid pairs: {len(pairs)}")

    # CONTROL CELLS
    ctrl_cells = adata.obs_names[adata.obs['gene'] == 'non-targeting']
    ctrl_counts = np.array((adata[ctrl_cells, :].X > expression_threshold).sum(axis=0)).flatten()
    ctrl_df = pd.DataFrame({'gene': adata.var_names, 'n_expr_ctrl': ctrl_counts})

    # TREATMENT CELLS
    trt_dfs = []
    for pert in tqdm([p for p in valid_perts if p != 'non-targeting']):
        pert_cells = adata.obs_names[adata.obs['gene'] == pert]
        counts = np.array((adata[pert_cells, :].X > expression_threshold).sum(axis=0)).flatten()
        trt_dfs.append(pd.DataFrame({
            'gene': adata.var_names,
            'perturbation': pert,
            'n_expr_trt': counts
        }))
    trt_counts_df = pd.concat(trt_dfs, ignore_index=True)

    # JOIN AND QC FLAG
    df = pairs.rename(columns={'effect_ensg': 'gene'})
    df = df.merge(trt_counts_df, on=['perturbation','gene'], how='left')
    df = df.merge(ctrl_df, on='gene', how='left')
    df = df.dropna(subset=['n_expr_trt','n_expr_ctrl'])
    df['QC_pass'] = np.where((df['n_expr_trt'] >= min_cells) & (df['n_expr_ctrl'] >= min_cells),
                             'Pass', 'Fail')
    return df

gene_QC = gene_qc(adata, pairs)

# -----------------------------
# MAP ENSG -> SYMBOL
# -----------------------------
mg = mygene.MyGeneInfo()
symbols = mg.querymany(gene_QC['gene'].unique(), scopes='ensembl.gene', fields='symbol', species='human')
symbol_map = {s['query']: s.get('symbol', None) for s in symbols}
gene_QC['effect'] = gene_QC['gene'].map(symbol_map)
# Don't filter perturbation == effect. Won't be able to estimate KD efficiency
# gene_QC = gene_QC[gene_QC['perturbation'] != gene_QC['effect']]

# -----------------------------
# SAVE RESULTS
# -----------------------------
gene_QC.to_csv(OUT_DIR / "gene_QC.csv", index=False)
print(f"QC done, saved to {OUT_DIR / 'gene_QC.csv'}")

