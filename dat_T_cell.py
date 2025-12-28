import scanpy as sc
import polars as pl
import pandas as pd
import sys
import numpy as np
from pathlib import Path


# Parse arguments
# h5ad_path = sys.argv[1]  # Path to h5ad file

h5ad_path = "/rds/project/rds-csoP2nj6Y6Y/biv22/data/GWCD4i.DE_stats.h5ad"

# Load data
adata = sc.read_h5ad(h5ad_path)

# Select only entries for gRNAs with successful on-target KD
on_target_adata = adata[adata.obs.ontarget_effect_category == "on-target KD"]

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
        "/rds/project/rds-csoP2nj6Y6Y/biv22/data/eqtl/cis_eQTLs_eQTLgen.txt",
        separator="\t",
    )
    .select(
        [
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
    pl.col("SNP"),
    pl.exclude("SNP").name.suffix(
        "_cis"
    ),  # append a _cis suffix to all columns except SNP
)
trans_eQTL = calculate_eQTL_beta(trans_eQTL, snp_frq_df).select(
    pl.col("SNP"),
    pl.exclude("SNP").name.suffix(
        "_trans"
    ),  # append a _trans suffix to all columns except SNP
)


merged_QTL = (
    cis_eQTL.join(
        trans_eQTL,
        on="SNP",
        how="inner",
    )
    .with_columns(y=(pl.col("FDR_trans") < 0.05).cast(pl.Int8))
    .rename({"GeneSymbol_cis": "cis_gene", "GeneSymbol_trans": "trans_gene"})
)

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

merged_QTL_lead = merged_QTL.join(
    lead_cis,
    on=["cis_gene", "SNP"],
    how="inner",
)

dat = merged_QTL_lead.join(
    DE_full,
    on=["cis_gene", "trans_gene"],
    how="inner",
)

dat.write_csv("/rds/project/rds-csoP2nj6Y6Y/biv22/data/pairs/full_dat_T_cell.csv")
