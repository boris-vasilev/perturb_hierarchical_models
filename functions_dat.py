import hail as hl
import polars as pl
import time


def harmonise_eqtl_1kg(eqtl_df: pl.DataFrame) -> pl.DataFrame:
    """
    Harmonise eQTL summary statistics with 1000G allele frequencies.

    Replicates R adjust_MAF + calculate_eQTL_beta logic:
        - join on SNP
        - align alleles
        - flip MAF if reversed
        - remove mismatches
        - compute SE and Beta
    """

    frq_1000g = (
        pl.scan_csv(
            "/rds/project/rds-csoP2nj6Y6Y/biv22/data/1000G_Phase3_frq/*.tsv",
            separator="\t",
        )
        .select(["SNP", "A1", "A2", "MAF"])
        .collect()
    )

    # join frequency data
    df = eqtl_df.lazy().join(frq_1000g.lazy(), on="SNP", how="inner")

    # allele alignment
    df = df.with_columns(
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
        .alias("MAF_aligned")
    )

    # drop mismatches
    df = df.filter(pl.col("MAF_aligned").is_not_null())

    # compute SE
    df = df.with_columns(
        (
            1.0
            / (
                2.0
                * pl.col("MAF_aligned")
                * (1.0 - pl.col("MAF_aligned"))
                * (pl.col("NrSamples") + pl.col("Zscore") ** 2)
            ).sqrt()
        ).alias("SE")
    )

    # compute Beta
    df = df.with_columns((pl.col("Zscore") * pl.col("SE")).alias("Beta"))

    return df.collect()


def harmonise_eqtl(eqtl_df: pl.DataFrame, AF: pl.DataFrame) -> pl.DataFrame:
    """
    Harmonise eQTL data with allele frequency data from eQTLgen. This includes:
        - Joining on SNP to get allele information
        - Determining if the alleles match, flip, or mismatch
        - Flipping Z-scores when alleles are flipped
        - Filtering out mismatches
    """
    # Read allele frequency data to calculate beta from Zscore and SE from NrSamples
    AF = (
        pl.scan_csv(
            "/rds/project/rds-csoP2nj6Y6Y/biv22/data/eqtl/AF_eQTLgen.txt",
            separator="\t",
            infer_schema_length=10000,
        )
        .select(
            [
                pl.col("SNP"),
                pl.col("AlleleA").alias("OtherAllele"),
                pl.col("AlleleB").alias("AssessedAllele"),
                pl.col("AlleleB_all").replace("NA", None).cast(pl.Float64).alias("AF"),
            ]
        )
        .collect()
    )
    df = eqtl_df.lazy().join(AF.lazy(), on="SNP", how="left", suffix="_AF")

    # allele relationship
    df = df.with_columns(
        pl.when(
            (pl.col("AssessedAllele") == pl.col("AssessedAllele_AF"))
            & (pl.col("OtherAllele") == pl.col("OtherAllele_AF"))
        )
        .then(pl.lit("match"))
        .when(
            (pl.col("AssessedAllele") == pl.col("OtherAllele_AF"))
            & (pl.col("OtherAllele") == pl.col("AssessedAllele_AF"))
        )
        .then(pl.lit("flip"))
        .otherwise(pl.lit("mismatch"))
        .alias("allele_status")
    )

    # flip Z-scores when alleles reversed
    df = df.with_columns(
        pl.when(pl.col("allele_status") == "flip")
        .then(-pl.col("Zscore"))
        .otherwise(pl.col("Zscore"))
        .alias("Zscore_aligned")
    )

    # optional but strongly recommended QC
    df = df.filter(pl.col("allele_status") != "mismatch")

    return df.collect()


def get_r2_between_snps(snp1, snp2):
    # Configure Hail for S3 access
    hl.init(
        spark_conf={
            "spark.jars.packages": "org.apache.hadoop:hadoop-aws:3.3.4",
            "spark.hadoop.fs.s3a.impl": "org.apache.hadoop.fs.s3a.S3AFileSystem",
            # Use ~/.aws/credentials
            "spark.hadoop.fs.s3a.aws.credentials.provider": "com.amazonaws.auth.DefaultAWSCredentialsProviderChain",
            # REQUIRED for pan-ukb
            "spark.hadoop.fs.s3a.requester.pays.enabled": "true",
        },
        idempotent=True,  # Run init only once
    )

    # Read variant indices table
    ht_idx = hl.read_table(
        "s3a://pan-ukb-us-east-1/ld_release/UKBB.EUR.ldadj.variant.ht"
    )

    # Load LD matrix
    bm = hl.linalg.BlockMatrix.read(
        "s3a://pan-ukb-us-east-1/ld_release/UKBB.EUR.ldadj.bm"
    )

    def flip_variant(v):
        chrom, pos, a1, a2 = v.split(":")
        return f"{chrom}:{pos}:{a2}:{a1}"

    flipped_snp1 = flip_variant(snp1)
    flipped_snp2 = flip_variant(snp2)

    # Variants provided as Chr:Pos:Ref:Alt strings which get parsed in the
    # appropriate locus alleles struct to reference the Hail table
    variants = [
        hl.parse_variant(snp1),
        hl.parse_variant(snp2),
        # To resolve issues with allele ordering (flipped ref/alt in LD matrix)
        # Because only one will exist in the LD matrix this will end up being
        # only 2 variants in the output
        hl.parse_variant(flipped_snp1),
        hl.parse_variant(flipped_snp2),
        # TODO: resolve ambiguities such as strand flips
        # TODO: drop those that cannot be resolved
    ]

    # Filter the variant index table to only the two variants of interest
    # and select their indices (to match in the LD BlockMatrix)
    ht_idx = ht_idx.filter(
        hl.literal(variants).contains(
            hl.struct(locus=ht_idx.locus, alleles=ht_idx.alleles)
        )
    )
    idx = ht_idx.idx.collect()

    # Get correlation from LD matrix.
    # Sort indices first because the LD matrix is symmetric but only the lower triangle is stored.
    i, j = sorted(idx)
    r = bm[i, j]

    return r**2


def compute_r2_for_pairs(
    df_pairs: pl.DataFrame,
    id_col: str = "id",
    lead_col: str = "id_lead",
):
    """
    df_pairs: Polars DF with columns [id, id_lead]
    returns: Polars DF with [id, id_lead, r2]
    """

    print("[LD] Initializing Hail")
    hl.init(
        spark_conf={
            "spark.driver.memory": "8g",
            "spark.executor.memory": "8g",
            "spark.executor.memoryOverhead": "2g",
            "spark.sql.shuffle.partitions": "200",
            "spark.default.parallelism": "200",
            "spark.jars.packages": "org.apache.hadoop:hadoop-aws:3.3.4",
            "spark.hadoop.fs.s3a.impl": "org.apache.hadoop.fs.s3a.S3AFileSystem",
            "spark.hadoop.fs.s3a.aws.credentials.provider": "com.amazonaws.auth.DefaultAWSCredentialsProviderChain",
            "spark.hadoop.fs.s3a.requester.pays.enabled": "true",
        },
        idempotent=True,
    )

    print("[LD] Loading variant index table")
    ht_idx = hl.read_table(
        "s3a://pan-ukb-us-east-1/ld_release/UKBB.EUR.ldadj.variant.ht"
    )

    print("[LD] Loading LD BlockMatrix")
    bm = hl.linalg.BlockMatrix.read(
        "s3a://pan-ukb-us-east-1/ld_release/UKBB.EUR.ldadj.bm"
    )

    def flip(v: str) -> str:
        c, p, a1, a2 = v.split(":")
        return f"{c}:{p}:{a2}:{a1}"

    print("[LD] Collecting pairs")
    pairs = df_pairs.to_dicts()

    all_variants = set()
    for r in pairs:
        all_variants.add(r[id_col])
        all_variants.add(r[lead_col])
        all_variants.add(flip(r[id_col]))
        all_variants.add(flip(r[lead_col]))

    print(f"[LD] Resolving {len(all_variants):,} variants to LD indices")

    t0 = time.time()

    ht_query = hl.Table.parallelize(
        [hl.parse_variant(v) for v in all_variants],
        key=["locus", "alleles"],
    )
    ht_query = ht_query.annotate(idx=ht_idx[ht_query.key].idx)

    rows = ht_query.collect()

    variant_to_idx = {}
    for r in rows:
        if r.idx is None:
            continue
        key = (
            f"{r.locus.contig}:"
            f"{r.locus.position}:"
            f"{r.alleles[0]}:"
            f"{r.alleles[1]}"
        )
        variant_to_idx[key] = int(r.idx)

    print(
        f"[LD] Resolved {len(variant_to_idx):,}/{len(all_variants):,} "
        f"variants in {time.time() - t0:.1f}s"
    )

    def get_idx(v: str):
        x = variant_to_idx.get(v)
        if x is not None:
            return x
        return variant_to_idx.get(flip(v))

    print("[LD] Computing LD")

    out = []
    total_pairs = len(pairs)

    for k, r in enumerate(pairs, start=1):
        i1 = get_idx(r[id_col])
        i2 = get_idx(r[lead_col])

        if i1 is None or i2 is None:
            r2 = None
        elif i1 == i2:
            r2 = 1.0
        else:
            i, j = sorted((i1, i2))
            r2 = float((bm[i, j]) ** 2)

        out.append({id_col: r[id_col], lead_col: r[lead_col], "r2": r2})

        if k % 1000 == 0:
            print(f"[LD] {k:,}/{total_pairs:,} pairs done")

    print("[LD] Finished LD computation")
    return pl.DataFrame(out)
