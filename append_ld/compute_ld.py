import os

# Increase Java Heap Size for PySpark
os.environ["PYSPARK_SUBMIT_ARGS"] = (
    "--driver-memory 6g "
    "--executor-memory 32g "
    "--conf spark.executor.memoryOverhead=10g "
    "pyspark-shell"
)
import sys
import hail as hl
import polars as pl
import time


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


input_file = sys.argv[1]

df = pl.read_csv(input_file)

out = compute_r2_for_pairs(df)

gene = input_file.split("/")[-1].replace(".csv", "")
out.write_csv(f"results/{gene}.ld.tsv")
