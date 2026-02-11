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

dat = merged_QTL.join(
    DE_full,
    on=["cis_gene", "trans_gene"],
    how="inner",
)

del DE_full, merged_QTL

dat = dat.with_columns(x=(pl.col("padj") < 0.05).cast(pl.Int8))

# Extract GCTA-COJO independent causal variants for each perturbed gene
cojo = (
    pl.read_csv(
        "/rds/project/rds-csoP2nj6Y6Y/biv22/data/eqtl/INTERVAL_eQTL_summary_statistics/cojo_independent_cis_signals.csv"
    )
    .select(
        pl.col("gene_name").alias("cis_gene"),
        pl.col("variant_id").alias("rsid_cojo"),
        pl.col("chr").cast(pl.Utf8).alias("chr"),
        pl.col("pos_b37").cast(pl.Int64).alias("pos_b37"),
        pl.col("effect_allele").alias("effect_allele"),
        pl.col("other_allele").alias("other_allele"),
    )
    # build the hail-friendly variant string
    .with_columns(
        pl.format(
            "{}:{}:{}:{}",
            pl.col("chr"),
            pl.col("pos_b37"),
            pl.col("other_allele"),
            pl.col("effect_allele"),
        ).alias("id_cojo")
    )
    .unique(subset=["cis_gene", "id_cojo"])
)

# 1) extract unique LD pairs
# Calculate r2 between SNP and COJO SNPs only for x1y1 pairs
# To reduce computation time and investigate only
df_ld_pairs = (
    dat.filter((pl.col("x") == 1) & (pl.col("y") == 1))
    .select(["cis_gene", "id"])
    .unique()
    .join(cojo.select(["cis_gene", "id_cojo"]), on="cis_gene", how="inner")
    .select(["id", "id_cojo"])
    .unique()
)

# 2) compute r2
r2_df = compute_r2_for_pairs(df_ld_pairs, lead_col="id_cojo")


r2_df.write_csv("/rds/project/rds-csoP2nj6Y6Y/biv22/data/pairs/gw_cojo_r2.csv")


# 3) get largest r2 for each gene across all COJO SNPs
r2_max = (
    dat.select(["cis_gene", "id"])
    .unique()
    .join(cojo.select(["cis_gene", "id_cojo"]), on="cis_gene", how="inner")
    .join(r2_df, on=["id", "id_cojo"], how="left")
    .group_by(["cis_gene", "id"])
    .agg(pl.max("r2").alias("r2_max_cojo"))
)
# get which COJO SNP achieved that max r2
r2_best = (
    dat.select(["cis_gene", "id"])
    .unique()
    .join(cojo.select(["cis_gene", "id_cojo"]), on="cis_gene", how="inner")
    .join(r2_df, on=["id", "id_cojo"], how="left")
    .join(r2_max, on=["cis_gene", "id"], how="inner")
    .filter(pl.col("r2") == pl.col("r2_max_cojo"))
    .group_by(["cis_gene", "id"])
    .first()  # if ties, just take the first one
    .select(["cis_gene", "id", "id_cojo", "r2_max_cojo"])
    .rename({"id_cojo": "id_cojo_best"})
)

# 4) join back to the full dataset
dat = dat.join(r2_best, on=["cis_gene", "id"], how="left")

dat.write_csv(
    "/rds/project/rds-csoP2nj6Y6Y/biv22/data/pairs/full_dat_GW_with_r2_cojo.csv"
)
