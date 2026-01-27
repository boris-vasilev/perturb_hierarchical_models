import hail as hl


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
