import hail as hl
import polars as pl
import time


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
              one row per unique LD query
    returns: Polars DF with [id, id_lead, r2]
    """

    # --------------------------------------------------
    # 1. Init Hail (memory is CRITICAL)
    # --------------------------------------------------
    print("[LD] Initializing Hail")
    hl.init(
        spark_conf={
            # memory
            "spark.driver.memory": "8g",
            "spark.executor.memory": "8g",
            "spark.executor.memoryOverhead": "2g",
            # avoid task explosion
            "spark.sql.shuffle.partitions": "200",
            "spark.default.parallelism": "200",
            # S3 / pan-UKB
            "spark.jars.packages": "org.apache.hadoop:hadoop-aws:3.3.4",
            "spark.hadoop.fs.s3a.impl": "org.apache.hadoop.fs.s3a.S3AFileSystem",
            "spark.hadoop.fs.s3a.aws.credentials.provider": "com.amazonaws.auth.DefaultAWSCredentialsProviderChain",
            "spark.hadoop.fs.s3a.requester.pays.enabled": "true",
        },
        idempotent=True,
    )

    # --------------------------------------------------
    # 2. Load LD resources
    # --------------------------------------------------
    print("[LD] Loading variant index table")
    ht_idx = hl.read_table(
        "s3a://pan-ukb-us-east-1/ld_release/UKBB.EUR.ldadj.variant.ht"
    )

    print("[LD] Loading LD BlockMatrix")
    bm = hl.linalg.BlockMatrix.read(
        "s3a://pan-ukb-us-east-1/ld_release/UKBB.EUR.ldadj.bm"
    )

    def flip(v):
        c, p, a1, a2 = v.split(":")
        return f"{c}:{p}:{a2}:{a1}"

    # --------------------------------------------------
    # 3. Collect unique variant strings
    # --------------------------------------------------
    print("[LD] Collecting pairs")
    pairs = df_pairs.to_dicts()

    all_variants = set()
    for r in pairs:
        all_variants |= {
            r[id_col],
            r[lead_col],
            flip(r[id_col]),
            flip(r[lead_col]),
        }

    print(f"[LD] Resolving {len(all_variants):,} variants to LD indices")

    # --------------------------------------------------
    # 4. Correct keyed lookup (NO join, NO scalar index)
    # --------------------------------------------------
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

    # --------------------------------------------------
    # 5. Group LD pairs by chromosome
    # --------------------------------------------------
    pairs_by_chr = {}
    for r in pairs:
        chrom = r[id_col].split(":")[0]
        pairs_by_chr.setdefault(chrom, []).append(r)

    print(f"[LD] Processing {len(pairs_by_chr)} chromosomes")

    # --------------------------------------------------
    # 6. Chromosome-wise LD computation (heap-safe)
    # --------------------------------------------------
    out = []
    total_done = 0
    total_pairs = len(pairs)

    for chrom, chrom_pairs in pairs_by_chr.items():
        print(f"[LD] chr{chrom}: {len(chrom_pairs):,} pairs")

        idxs = set()
        pair_indices = []

        for r in chrom_pairs:
            i1 = variant_to_idx.get(r[id_col]) or variant_to_idx.get(flip(r[id_col]))
            i2 = variant_to_idx.get(r[lead_col]) or variant_to_idx.get(
                flip(r[lead_col])
            )
            pair_indices.append((i1, i2))
            if i1 is not None:
                idxs.add(i1)
            if i2 is not None:
                idxs.add(i2)

        if not idxs:
            continue

        idxs = sorted(idxs)

        print(f"[LD] chr{chrom}: extracting {len(idxs)}×{len(idxs)} LD block")
        t0 = time.time()

        bm_sub = bm.filter_rows(idxs).filter_cols(idxs).to_numpy()
        idx_to_pos = {idx: i for i, idx in enumerate(idxs)}

        print(f"[LD] chr{chrom}: LD loaded in {time.time() - t0:.1f}s")

        for (i1, i2), r in zip(pair_indices, chrom_pairs):
            if i1 is None or i2 is None:
                r2 = None
            elif i1 == i2:
                r2 = 1.0
            else:
                r2 = bm_sub[idx_to_pos[i1], idx_to_pos[i2]] ** 2

            out.append(
                {
                    id_col: r[id_col],
                    lead_col: r[lead_col],
                    "r2": float(r2) if r2 is not None else None,
                }
            )

        total_done += len(chrom_pairs)
        print(
            f"[LD] Progress: {total_done}/{total_pairs} "
            f"({100 * total_done / total_pairs:.1f}%)"
        )

    print("[LD] Finished LD computation")

    return pl.DataFrame(out)
