import os

# Increase Java Heap Size for PySpark
os.environ["PYSPARK_SUBMIT_ARGS"] = (
    "--driver-memory 6g "
    "--executor-memory 32g "
    "--conf spark.executor.memoryOverhead=10g "
    "pyspark-shell"
)


import scanpy as sc
import polars as pl
import pandas as pd
import numpy as np
from pathlib import Path
from functions_LD import compute_r2_for_pairs

# Parse arguments
# h5ad_path = sys.argv[1]  # Path to h5ad file

h5ad_path = "/rds/project/rds-csoP2nj6Y6Y/biv22/data/GWCD4i.DE_stats.h5ad"

# Load data
adata = sc.read_h5ad(h5ad_path)

# Select only entries for gRNAs with successful on-target KD
on_target_adata = adata[adata.obs.ontarget_effect_category == "on-target KD"]

del adata

perturbed_genes = on_target_adata.obs.target_contrast_gene_name.unique()
expressed_genes = on_target_adata.var.gene_name.unique()

obs_df = pl.from_pandas(
    on_target_adata.obs[
        [
            "target_contrast",
            "target_contrast_gene_name",
            "culture_condition",
            "ontarget_effect_size",
            "target_baseMean",
        ]
    ].reset_index(names="obs_id")
).with_columns(pl.col("target_contrast_gene_name").cast(pl.String))

var_df = pl.from_pandas(
    on_target_adata.var[
        [
            "gene_ids",
            "gene_name",
        ]
    ].reset_index(names="var_id")
)


def layer_to_df(adata, layer_name, value_name):
    mat = adata.layers[layer_name]

    df = pl.DataFrame(
        {
            "obs_id": np.repeat(adata.obs_names, mat.shape[1]),
            "var_id": np.tile(adata.var_names, mat.shape[0]),
            value_name: mat.ravel(order="C"),
        }
    )
    return df


adj_p = layer_to_df(on_target_adata, "adj_p_value", "padj")
baseMean = layer_to_df(on_target_adata, "baseMean", "base_DESeq2")
lfcSE = layer_to_df(on_target_adata, "lfcSE", "lfcSE")
log_fc = layer_to_df(on_target_adata, "log_fc", "logFC")

DE_df = pl.concat(
    [
        adj_p,
        baseMean.select("base_DESeq2"),
        lfcSE.select("lfcSE"),
        log_fc.select("logFC"),
    ],
    how="horizontal",
)

## Concatenate all the DE summary statistics in a single df.
## Because order is preserved from the anndata, we don't need to use joins. Only concatenate/horizontal stacking
n_obs = obs_df.height  # number of perturbs
n_var = var_df.height  # number of expressed genes

obs_idx = pl.Series(
    "obs_idx",
    [i for i in range(n_obs) for _ in range(n_var)],
)

var_idx = pl.Series(
    "var_idx",
    list(range(n_var)) * n_obs,
)

DE_full = DE_df.with_columns(
    [
        obs_df["culture_condition"].gather(obs_idx).alias("culture_condition"),
        obs_df["target_contrast"].gather(obs_idx).alias("cis_ensg"),
        obs_df["target_contrast_gene_name"].gather(obs_idx).alias("cis_gene"),
        obs_df["ontarget_effect_size"].gather(obs_idx).alias("perturb_eff"),
        obs_df["target_baseMean"].gather(obs_idx).alias("cis_baseMean"),
        var_df["gene_ids"].gather(var_idx).alias("trans_ensg"),
        var_df["gene_name"].gather(var_idx).alias("trans_gene"),
    ]
)

# Read cis-eQTLs
cis_eQTL = (
    pl.scan_csv(
        "/rds/project/rds-csoP2nj6Y6Y/biv22/data/eqtl/cis_eQTLs_eQTLgen.tsv",
        separator="\t",
    )
    .select(
        [
            pl.format(
                "{}:{}:{}:{}",
                pl.col("SNPChr"),
                pl.col("SNPPos"),
                pl.col("OtherAllele"),
                pl.col("AssessedAllele"),
            ).alias("id"),
            "SNP",
            "Pvalue",
            "FDR",
            "GeneSymbol",
            "Zscore",
            "NrSamples",
            "AssessedAllele",
            "OtherAllele",
        ],
    )
    .filter(
        pl.col("GeneSymbol").is_in(
            perturbed_genes
        ),  # Only cis-eQTLs of perturbed genes
        pl.col("FDR") < 0.05,  # Only significant cis-eQTLs
    )
    .collect()
)
# Read trans-eQTLs
trans_eQTL = (
    pl.scan_csv(
        "/rds/project/rds-csoP2nj6Y6Y/biv22/data/eqtl/trans_eQTLs_eQTLgen.txt",
        separator="\t",
    )
    .select(
        [
            pl.format(
                "{}:{}:{}:{}",
                pl.col("SNPChr"),
                pl.col("SNPPos"),
                pl.col("OtherAllele"),
                pl.col("AssessedAllele"),
            ).alias("id"),
            "SNP",
            "Pvalue",
            "FDR",
            "GeneSymbol",
            "Zscore",
            "NrSamples",
            "AssessedAllele",
            "OtherAllele",
        ]
    )
    .filter(
        pl.col("GeneSymbol").is_in(
            expressed_genes
        )  # Only trans-eQTLs of expressed genes
    )
    .collect()
)

frq_dir = Path("/rds/project/rds-csoP2nj6Y6Y/biv22/data/1000G_Phase3_frq")

dfs = []
for f in frq_dir.glob("*.frq"):
    pdf = pd.read_fwf(f)  # fixed-width parser (fread-like)
    dfs.append(pl.from_pandas(pdf))

snp_frq_df = pl.concat(dfs, how="vertical")


def adjust_MAF(eqtl_df: pl.DataFrame, snp_frq_df: pl.DataFrame) -> pl.DataFrame:
    return (
        eqtl_df.join(snp_frq_df, on="SNP", how="inner")
        .with_columns(
            pl.when(
                (pl.col("AssessedAllele") == pl.col("A1"))
                & (pl.col("OtherAllele") == pl.col("A2"))
            )
            .then(pl.col("MAF"))
            .when(
                (pl.col("AssessedAllele") == pl.col("A2"))
                & (pl.col("OtherAllele") == pl.col("A1"))
            )
            .then(1.0 - pl.col("MAF"))
            .otherwise(None)
            .alias("MAF")
        )
        .filter(pl.col("MAF").is_not_null())
    )


def calculate_eQTL_beta(
    eqtl_df: pl.DataFrame, snp_frq_df: pl.DataFrame
) -> pl.DataFrame:
    eqtl_df = adjust_MAF(eqtl_df, snp_frq_df)

    return (
        eqtl_df.with_columns(
            SE=1.0
            / (
                2.0
                * pl.col("MAF")
                * (1.0 - pl.col("MAF"))
                * (pl.col("NrSamples") + pl.col("Zscore") ** 2)
            ).sqrt(),
        )
        .with_columns(Beta=pl.col("Zscore") * pl.col("SE"))
        .select(
            [
                "id",
                "SNP",
                "GeneSymbol",
                "Pvalue",
                "FDR",
                "Beta",
                "SE",
                "Zscore",
                "MAF",
            ]
        )
    )


cis_eQTL = calculate_eQTL_beta(cis_eQTL, snp_frq_df).select(
    pl.col("id"),
    pl.col("SNP"),
    pl.exclude("SNP").name.suffix(
        "_cis"
    ),  # append a _cis suffix to all columns except SNP
)

# Identify lead SNPs (smallest p-value)
lead_snps = (
    cis_eQTL.sort(
        ["GeneSymbol_cis", "Pvalue_cis", pl.col("Beta_cis").abs()],
        descending=[False, False, True],
    )
    .group_by("GeneSymbol_cis")
    .first()
    .select(
        pl.col("GeneSymbol_cis").alias("cis_gene"),
        pl.col("SNP").alias("rsid_lead"),
        pl.col("id").alias("id_lead"),
    )
)

trans_eQTL = calculate_eQTL_beta(trans_eQTL, snp_frq_df).select(
    pl.col("SNP"),
    pl.exclude("SNP").name.suffix(
        "_trans"
    ),  # append a _trans suffix to all columns except SNP
)

del dfs, snp_frq_df

merged_QTL = (
    cis_eQTL.join(
        trans_eQTL,
        on="SNP",
        how="inner",
    )
    .with_columns(y=(pl.col("FDR_trans") < 0.05).cast(pl.Int8))
    .rename({"GeneSymbol_cis": "cis_gene", "GeneSymbol_trans": "trans_gene"})
)

del cis_eQTL, trans_eQTL

lead_cis = (
    merged_QTL.select(["cis_gene", "SNP", "Pvalue_cis", "Beta_cis"])
    .sort(
        ["cis_gene", "Pvalue_cis", pl.col("Beta_cis").abs()],
        descending=[False, False, True],
    )
    .group_by("cis_gene")
    .first()
    .select(["cis_gene", "SNP"])
)

merged_QTL = merged_QTL.join(
    lead_cis,
    on=["cis_gene", "SNP"],
    how="inner",
)

dat = merged_QTL.join(
    DE_full,
    on=["cis_gene", "trans_gene"],
    how="inner",
)

del DE_full, merged_QTL

dat = dat.with_columns(x=(pl.col("padj") < 0.05).cast(pl.Int8))

# Join true lead SNP information to estimate r2 with true lead SNPs
dat = dat.join(
    lead_snps,
    on="cis_gene",
    how="inner",
)

# 1) extract unique LD pairs
# Calculate r2 between SNP and lead SNP only for x1y1 pairs
# To reduce computation time and investigate only
df_ld_pairs = (
    dat.filter(pl.col("x") == 1, pl.col("y") == 1).select(["id", "id_lead"]).unique()
)

# 2) compute r2
r2_df = compute_r2_for_pairs(df_ld_pairs)

# 3) join back to the full dataset
dat = dat.join(
    r2_df,
    on=["id", "id_lead"],
    how="left",
)

dat_split = dat.partition_by("culture_condition", as_dict=True)

for condition, condition_dat in dat_split.items():
    condition_dat.write_csv(
        f"/rds/project/rds-csoP2nj6Y6Y/biv22/data/pairs/full_dat_T_cell_{condition[0]}_with_r2.csv"
    )
