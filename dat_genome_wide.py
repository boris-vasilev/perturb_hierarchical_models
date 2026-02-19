import os
import argparse
import polars as pl

# Increase Java Heap Size for PySpark
os.environ["PYSPARK_SUBMIT_ARGS"] = (
    "--driver-memory 6g "
    "--executor-memory 32g "
    "--conf spark.executor.memoryOverhead=10g "
    "pyspark-shell"
)

from pathlib import Path
from functions_dat import compute_r2_for_pairs, harmonise_eqtl, harmonise_eqtl_1kg

parser = argparse.ArgumentParser()
parser.add_argument("--screen", type=str, required=True, help="Perturb-seq screen")
parser.add_argument("--output", type=str, required=True, help="Output file name")
parser.add_argument("--ld", action="store_true")
parser.add_argument("--no-ld", dest="ld", action="store_false")
parser.set_defaults(ld=True)
parser.add_argument("--beta", action="store_true")
parser.add_argument("--no-beta", dest="beta", action="store_false")
parser.set_defaults(beta=True)
parser.add_argument(
    "--af",
    type=str,
    help="Allele frequency source",
    choices=["eQTLgen", "1KG"],
    default="eQTLgen",
)
args = parser.parse_args()
screen = args.screen

eqtl_study = "eQTLgen"

DE_output_dir = "/rds/project/rds-csoP2nj6Y6Y/biv22/data/perturb"
QC_dir = "/rds/project/rds-csoP2nj6Y6Y/biv22/data/perturb_QC"
eqtl_dir = "/rds/project/rds-csoP2nj6Y6Y/biv22/data/eqtl"
# Input files
QC_file = f"{QC_dir}/{screen}/gene_QC.csv"
base_expression_file = f"{QC_dir}/{screen}/base_expression.csv"
# Output file
output_file = f"/rds/project/rds-csoP2nj6Y6Y/biv22/data/pairs/{args.output}"


### Read Perturb-seq data
# Read QC file to find pairs passing pairwise QC
QC = (
    pl.scan_csv(QC_file)
    .filter(pl.col("n_expr_trt") > 3, pl.col("n_expr_ctrl") > 3)
    .select(pl.col("perturbation"), pl.col("gene"), pl.col("effect"))
    .collect()
)

print("[1/9] Read DE stats")

DE_stats = (
    pl.scan_csv(
        f"{DE_output_dir}/{screen}/*.tsv",
        separator="\t",
        include_file_paths="filepath",
        schema_overrides={
            "gene": pl.Utf8,
            "log2FoldChange": pl.Float64,
            "lfcSE": pl.Float64,
            "padj": pl.Float64,
        },
        null_values=["NA"],
    )
    .select(
        [
            pl.col("filepath").str.extract(r"([^/]+)\.tsv$").alias("cis_gene"),
            pl.col("gene").alias("trans_gene"),
            pl.col("log2FoldChange").alias("logFC"),
            pl.col("lfcSE"),
            pl.col("padj"),
            pl.when(pl.col("padj") < 0.05)
            .then(1)
            .otherwise(0)
            .cast(pl.Int8)
            .alias("x"),
        ]
    )
    .collect()
)

# Filter based on QC (pairwise)
DE_stats = DE_stats.join(
    QC,
    left_on=["cis_gene", "trans_gene"],
    right_on=["perturbation", "gene"],
    how="inner",
).select(
    pl.col("cis_gene"),
    pl.col("trans_gene").alias("trans_gene_ensg"),
    pl.col("logFC"),
    pl.col("lfcSE"),
    pl.col("padj"),
    pl.col("x"),
    pl.col("effect").alias("trans_gene"),
)
## Add perturbation efficiency and efficiency SE
# get cis (self) perturbation effects
cis_effects = (
    DE_stats.lazy()
    .filter(pl.col("cis_gene") == pl.col("trans_gene"))
    .select(
        [
            pl.col("cis_gene"),
            pl.col("logFC").alias("perturb_eff"),
            pl.col("lfcSE").alias("perturb_eff_se"),
        ]
    )
)

# join back onto full table
DE_stats = DE_stats.lazy().join(cis_effects, on="cis_gene", how="left").collect()

# Join base expression data (from the expression in NT cells)
base_expression = pl.read_csv(base_expression_file).select(
    [
        pl.col("gene_id").alias("trans_gene_ensg"),
        # pl.col("base_mean_per_cell"),
        # pl.col("log1p_base_mean_per_cell"),
        # pl.col("base_cpm"),
        # pl.col("log1p_base_cpm"),
        pl.col("base_deseq2"),
        pl.col("log1p_base_deseq2"),
    ]
)
DE_stats = DE_stats.join(base_expression, on="trans_gene_ensg", how="left")


### Read eQTL data
print("[2/9] Read cis-eQTLs")
# Get list of cis/trans genes from the DE stats to filter eQTLs and not load the whole summary stats
cis_genes = DE_stats.select("cis_gene").unique().to_series().to_list()
trans_genes = DE_stats.select("trans_gene").unique().to_series().to_list()
trans_genes_ensg = DE_stats.select("trans_gene_ensg").unique().to_series().to_list()

# Read cis-eQTLs
cis_eQTL = (
    pl.scan_csv(
        f"{eqtl_dir}/cis_eQTLs_{eqtl_study}.tsv",
        separator="\t",
    )
    # .filter(
    #     ~(
    #         (pl.col("SNPChr") == 6)
    #         & (pl.col("SNPPos").is_between(28_477_797, 33_448_354))
    #     )  # Exclude MHC region
    # )
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
            "Gene",
            "GeneSymbol",
            "Zscore",
            "NrSamples",
            "AssessedAllele",
            "OtherAllele",
        ],
    )
    .filter(
        pl.col("GeneSymbol").is_in(cis_genes),  # Only cis-eQTLs of perturbed genes
        pl.col("FDR") < 0.05,  # Only significant cis-eQTLs
        # Filter out the IFNAR2-IL10RB readthrough gene (ENSG00000249624) which has the same symbol as IFNAR2.
        # Leaving it in introduces duplicates for IFNAR2.
        # Removing it from cis-eQTLs is consistent with what the Perturb-seq is -- perturbation of IFNAR2 not the readthrough transcript.
        # pl.col("Gene") != "ENSG00000249624",
    )
    .collect()
)

print("[3/9] Read trans-eQTLs")
# Read trans-eQTLs
trans_eQTL = (
    pl.scan_csv(
        f"{eqtl_dir}/trans_eQTLs_{eqtl_study}.txt",
        separator="\t",
    )
    .select(
        [
            # pl.format(
            #     "{}:{}:{}:{}",
            #     pl.col("SNPChr"),
            #     pl.col("SNPPos"),
            #     pl.col("OtherAllele"),
            #     pl.col("AssessedAllele"),
            # ).alias("id"),
            "SNP",
            "Pvalue",
            "FDR",
            "Gene",
            "GeneSymbol",
            "Zscore",
            "NrSamples",
            "AssessedAllele",
            "OtherAllele",
        ]
    )
    .filter(
        pl.col("SNP").is_in(cis_eQTL["SNP"]),
        pl.col("Gene").is_in(trans_genes_ensg),  # Only trans-eQTLs of expressed genes
    )
    .collect()
)

if args.beta:
    print("[4/9] Calculate Beta and SE for eQTLs")
    # Harmonise eQTLs with allele frequency data

    if args.af == "eQTLgen":
        cis_eQTL = harmonise_eqtl(cis_eQTL)
        trans_eQTL = harmonise_eqtl(trans_eQTL)

        # Calculate Beta and SE from Z-score and NrSamples
        cis_eQTL = cis_eQTL.with_columns(
            SE=1.0
            / (
                2.0
                * pl.col("AF")
                * (1.0 - pl.col("AF"))
                * (pl.col("NrSamples") + pl.col("Zscore_aligned") ** 2)
            ).sqrt(),
        ).with_columns(Beta=pl.col("Zscore_aligned") * pl.col("SE"))

        trans_eQTL = trans_eQTL.with_columns(
            SE=1.0
            / (
                2.0
                * pl.col("AF")
                * (1.0 - pl.col("AF"))
                * (pl.col("NrSamples") + pl.col("Zscore_aligned") ** 2)
            ).sqrt(),
        ).with_columns(Beta=pl.col("Zscore_aligned") * pl.col("SE"))
    elif args.af == "1KG":
        cis_eQTL = harmonise_eqtl_1kg(cis_eQTL)
        trans_eQTL = harmonise_eqtl_1kg(trans_eQTL)

# # Add cis/trans suffix to columns (except SNP) to avoid name clashes when merging
cis_eQTL = cis_eQTL.select(
    pl.col("id"),
    pl.col("SNP"),
    pl.exclude(["id", "SNP"]).name.suffix(
        "_cis"
    ),  # append a _cis suffix to all columns except SNP and id
)

# # Identify lead SNPs (smallest p-value) (before merge with trans-eQTLs. Select the lead SNP!)
lead_snps = (
    cis_eQTL.sort(
        [
            "GeneSymbol_cis",
            "Pvalue_cis",
            pl.col("Beta_cis" if args.beta else "Zscore_cis").abs(),
        ],
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

trans_eQTL = trans_eQTL.select(
    pl.col("SNP"),
    pl.exclude("SNP").name.suffix(
        "_trans"
    ),  # append a _trans suffix to all columns except SNP
)
print("[5/9] Merge cis and trans eQTLs to form eQTL pairs")
merged_QTL = (
    cis_eQTL.join(
        trans_eQTL,
        on="SNP",
        how="inner",
    )
    .with_columns(y=(pl.col("FDR_trans") < 0.05).cast(pl.Int8))
    .rename(
        {
            "GeneSymbol_cis": "cis_gene",
            "GeneSymbol_trans": "trans_gene",
            "Gene_trans": "trans_gene_ensg",
        }
    )
)
del cis_eQTL, trans_eQTL

print(
    "[6/9] Select only the most significant cis-eSNP tested for trans effects (WARNING: might not be lead SNP.)"
)
# Select top SNPs per cis gene (smallest p-value) - this is the most significant
# cis-eSNP tested for trans effects, but not necessarily the lead cis-eSNP (smallest p-value overall) as that one might not have been tested for trans effects.
# To get the lead cis-eSNP, we would need to go back to the full cis-eQTL summary stats, find the lead SNP per gene, and check if it was tested for trans effects
# NOTE: For eQTLgen -- the lead SNP is often not trans-tested. For INTERVAL -- the lead SNPs are the ones tested
# selected_snps = (
#     merged_QTL.with_columns(
#         abs_cis_eff=pl.col("Beta_cis" if args.beta else "Zscore_cis").abs()
#     )
#     .sort(["cis_gene", "Pvalue_cis", "abs_cis_eff"], descending=[False, False, True])
#     .unique(subset=["cis_gene", "SNP"], keep="first")
#     .unique(subset=["cis_gene"], keep="first")
#     .select(["cis_gene", "SNP"])
# )

selected_snps = (
    merged_QTL.group_by("cis_gene")
    .agg([pl.min("Pvalue_cis").alias("min_p")])
    .join(merged_QTL, on="cis_gene")
    .filter(pl.col("Pvalue_cis") == pl.col("min_p"))
    .with_columns(abs_cis_eff=pl.col("Beta_cis" if args.beta else "Zscore_cis").abs())
    .sort(["cis_gene", "abs_cis_eff"], descending=[False, True])
    .group_by("cis_gene")
    .first()
    .select(["cis_gene", "SNP"])
)


# Merge back onto full table to keep all columns but only the top SNP per cis gene
merged_QTL = merged_QTL.join(selected_snps, on=["cis_gene", "SNP"], how="inner")

print("[7/9] Merge eQTL pairs with perturbation pairs")
dat = merged_QTL.join(
    DE_stats,
    on=["cis_gene", "trans_gene_ensg"],
    how="inner",
)

del DE_stats, merged_QTL

if args.ld:
    # Calculate r2 between cis-gene lead SNP (from)
    print(
        "[8/9] Calculate r2 between selected top cis-eSNP and lead cis-eSNP from eQTLgen (WARNING: for x1y1 pairs only)"
    )

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
        dat.filter(pl.col("x") == 1, pl.col("y") == 1)
        .select(["id", "id_lead"])
        .unique()
    )

    # 2) compute r2
    r2_df_lead = compute_r2_for_pairs(df_ld_pairs).rename({"r2": "r2_lead"})

    # 3) join back to the full dataset
    dat = dat.join(
        r2_df_lead,
        on=["id", "id_lead"],
        how="left",
    )

    # Calculate r2 for each cis-gene selected SNP (top trans tested from eQTLgen)
    # with the cis-gene independent cis signals from INTERVAL (GCTA-COJO)
    print(
        "[9/9] Calculate r2 between selected top cis-eSNP and independent cis signals from INTERVAL (GCTA-COJO) (WARNING: for x1y1 pairs only)"
    )
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
    r2_df_cojo = compute_r2_for_pairs(df_ld_pairs, lead_col="id_cojo").rename(
        {"r2": "r2_cojo"}
    )

    # 3) get largest r2 for each gene across all COJO SNPs
    r2_max = (
        dat.select(["cis_gene", "id"])
        .unique()
        .join(cojo.select(["cis_gene", "id_cojo"]), on="cis_gene", how="inner")
        .join(r2_df_cojo, on=["id", "id_cojo"], how="left")
        .group_by(["cis_gene", "id"])
        .agg(
            pl.max("r2_cojo").alias("r2_max_cojo"),
            pl.count("id_cojo").alias("n_cojo_`signals"),
        )
    )
    # get which COJO SNP achieved that max r2
    r2_best = (
        dat.select(["cis_gene", "id"])
        .unique()
        .join(cojo.select(["cis_gene", "id_cojo"]), on="cis_gene", how="inner")
        .join(r2_df_cojo, on=["id", "id_cojo"], how="left")
        .join(r2_max, on=["cis_gene", "id"], how="inner")
        .filter(pl.col("r2_cojo") == pl.col("r2_max_cojo"))
        .group_by(["cis_gene", "id"])
        .first()  # if ties, just take the first one
        .select(["cis_gene", "id", "id_cojo", "r2_max_cojo"])
        .rename({"id_cojo": "id_cojo_best"})
    )

    # 4) join back to the full dataset
    dat = dat.join(r2_best, on=["cis_gene", "id"], how="left")

print("✔ Writing final output...")

dat.write_csv(output_file)
