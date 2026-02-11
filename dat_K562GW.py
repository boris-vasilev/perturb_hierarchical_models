import os

# Increase Java Heap Size for PySpark
os.environ["PYSPARK_SUBMIT_ARGS"] = (
    "--driver-memory 6g "
    "--executor-memory 32g "
    "--conf spark.executor.memoryOverhead=10g "
    "pyspark-shell"
)


import polars as pl
import pandas as pd
from pathlib import Path
from functions_LD import compute_r2_for_pairs


DE_full = pl.read_csv(
    "/rds/project/rds-csoP2nj6Y6Y/biv22/data/pairs/full_dat_GW.csv",
    infer_schema_length=100000,
).select(
    pl.col("perturb").alias("cis_gene"),
    pl.col("effect").alias("trans_gene"),
    pl.col("perturb_eff"),
    pl.col("perturb_eff_se"),
    pl.col("logFC"),
    pl.col("lfcSE"),
    pl.col("padj"),
    pl.col("x"),
    pl.col("base_DESeq2"),
    pl.col("base_log1p_DESeq2"),
)


perturbed_genes = DE_full.select("cis_gene").unique().to_series().to_list()
expressed_genes = DE_full.select("trans_gene").unique().to_series().to_list()

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

dat.write_csv("/rds/project/rds-csoP2nj6Y6Y/biv22/data/pairs/full_dat_GW_with_r2.csv")
