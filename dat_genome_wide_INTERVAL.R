library(tidyverse)
library(data.table)
library(parallel)
library(ashr)
library(argparse)
library(org.Hs.eg.db)
source("functions_create_perturb_df.R")

parser <- ArgumentParser()
parser$add_argument("--cores", type = "numeric", help = "Number of cores", required = TRUE)

args <- parser$parse_args()

n_cores <- args$cores

cells <- c("K562_GenomeWide")

data_dir <- "/rds/project/rds-csoP2nj6Y6Y/biv22/data"
DE_dir <- file.path(data_dir, "perturb")
QC_dir <- file.path(data_dir, "perturb_QC")
pairs_dir <- file.path(data_dir, "pairs")

output_file <- file.path(pairs_dir, "full_dat_GW_INERVAL.csv")

cell_DE_dirs <- setNames(file.path(DE_dir, cells), cells)

cell_DEG_files <- lapply(cell_DE_dirs, function(dir) {list.files(dir, full.names = TRUE)})

gene_QC_files <- setNames(file.path(QC_dir, cells, "gene_QC.csv"), cells)
base_expression_files <- setNames(file.path(QC_dir, cells, "base_expression.csv"), cells)

gene_QC <- mclapply(gene_QC_files, fread, mc.cores = length(cells))
base_expressions <- mclapply(base_expression_files, fread, mc.cores = length(cells))

# Select perturbation-gene pairs that pass QC (> 7 expressing cells in perturb and control groups)
QC_pairs_pass <- mclapply(gene_QC, function(QC) {
  QC %>%
    filter(n_expr_trt > 3,  n_expr_ctrl > 3) %>%
    dplyr::select(perturbation, gene, effect)
}, mc.cores = length(cells))


print("Extracting DEGs")

# number of cores per cell
cores_per_cell <- max(1, floor(n_cores / length(cells)))

results_per_cell <- mclapply(seq_along(cell_DEG_files), function(i) {
  cell <- cells[i]
  DEG_files <- cell_DEG_files[[i]]
  message("Processing cell: ", cell, " with ", cores_per_cell, " cores")

  QC_pass <- QC_pairs_pass[[i]]
  base_exp <- base_expressions[[i]]

  perturb_effects <- mclapply(DEG_files, function(file) {
    tryCatch({
      perturbed_gene <- basename(file) %>% sub("\\.tsv$", "", .)
      message("Reading file: ", file)

      DEGs <- fread(file)
      required_cols <- c("gene", "log2FoldChange", "lfcSE", "padj")
      if (!all(required_cols %in% colnames(DEGs))) {
        warning("Skipping ", file, " (missing columns)")
        return(NULL)
      }

      DEGs <- DEGs %>%
        inner_join(QC_pass %>% filter(perturbation == perturbed_gene), by = "gene")
      if (nrow(DEGs) == 0) return(NULL)

      ash_fit <- tryCatch(
        ashr::ash(beta = DEGs$log2FoldChange, se = DEGs$lfcSE, method = "shrink"),
        error = function(e) NULL
      )
      if (is.null(ash_fit)) return(NULL)

      # Add baseline expression estimates (mean UMI/cell, DESeq2 median-of-ratios, and CPM, + log1p transformed versions)
      DEGs <- merge(DEGs, base_exp, by.x = "gene", by.y = "gene_id", all.x = TRUE)

      pert_eff_value <- DEGs$log2FoldChange[DEGs$effect == perturbed_gene]
      pert_eff_se <- DEGs$lfcSE[DEGs$effect == perturbed_gene]
      # pert_eff_value <- 1 - 2^(DEGs$log2FoldChange[DEGs$effect == perturbed_gene])
      # if (length(pert_eff_value) == 0) pert_eff_value <- 1
      DEGs$perturb_eff <- pert_eff_value
      DEGs$perturb_eff_se <- pert_eff_se

      data.frame(
        perturb = perturbed_gene,
        gene = DEGs$gene,
        effect = DEGs$effect,
        # base_mean_per_cell = DEGs$base_mean_per_cell,
        # base_log1p_mean_per_cell = DEGs$log1p_base_mean_per_cell,
        # base_CPM = DEGs$base_cpm,
        # base_log1p_CPM = DEGs$log1p_base_cpm,
        base_DESeq2 = DEGs$base_deseq2,
        base_log1p_DESeq2 = DEGs$log1p_base_deseq2,
        perturb_eff = DEGs$perturb_eff,
        perturb_eff_se = DEGs$perturb_eff_se,
        logFC = DEGs$log2FoldChange,
        lfcSE = DEGs$lfcSE,
        padj = DEGs$padj,
        x = ifelse(DEGs$padj < 0.05, 1, 0),
        ash_logFC = ash_fit$result$PosteriorMean,
        ash_lfcSE = ash_fit$result$PosteriorSD,
        ash_svalue = ash_fit$result$svalue,
        # ash_x = ifelse(ash_fit$result$svalue < 0.05, 1, 0),
        screen = cell
      )
    }, error = function(e) {
      warning("Error in file ", file, ": ", conditionMessage(e))
      return(NULL)
    })
  }, mc.cores = cores_per_cell, mc.preschedule = FALSE)

  perturb_effects <- Filter(Negate(is.null), perturb_effects)  # remove NULLs
  if (length(perturb_effects) == 0) return(NULL)

  bind_rows(perturb_effects)
}, mc.cores = length(cells))

names(results_per_cell) <- cells

perturb_summary_stats <- rbindlist(results_per_cell)

####### PROCESS EQTLS
### -- Load and filter cis-eQTLs -- ###
message("[1/7] Reading cis-eQTLs")
cis_eQTL <- fread("/rds/project/rds-csoP2nj6Y6Y/biv22/data/eqtl/cis_eQTLs_INTERVAL.tsv",
                  sep = "\t",
                  # select = c("SNP", "Pvalue", "FDR", "GeneSymbol", "Zscore", "NrSamples", "AssessedAllele", "OtherAllele")
                  select = c("variant_rsid","phenotype_id", "gene_name", "pval_nominal", "qval", "slope", "slope_se")
                  ) %>%
                  dplyr::rename(SNP = variant_rsid,
                                ENSG = phenotype_id,
                                GeneSymbol = gene_name,
                                Pvalue = pval_nominal,
                                FDR = qval,
                                Beta = slope,
                                SE = slope_se)
message("  ✔ cis-eQTLs loaded: ", nrow(cis_eQTL))

message("[2/7] Selecting significant cis-eSNPs and calculating beta")
setkey(cis_eQTL, GeneSymbol)
cis_eSNPs <- cis_eQTL %>%
  as.data.table %>%
  filter(GeneSymbol %in% perturb_summary_stats$perturb)

message("  ✔ Significant perturbation cis-eSNPs found: ", nrow(cis_eSNPs))

message("[3/7] Reading trans-eQTLs")
trans_eQTL <- fread("/rds/project/rds-csoP2nj6Y6Y/biv22/data/eqtl/trans_eQTLs_INTERVAL.tsv",
                    sep = "\t",
                    # select = c("SNP", "Pvalue", "FDR", "GeneSymbol", "Zscore", "NrSamples", "AssessedAllele", "OtherAllele"),
                    select = c("variant_id","phenotype_id", "gene_name", "pval", "b", "b_se")) %>%
                    dplyr::rename(SNP = variant_id,
                                  ENSG = phenotype_id,
                                  Pvalue = pval,
                                  Beta = b,
                                  SE = b_se)
message("  ✔ trans-eQTLs loaded: ", nrow(trans_eQTL))

trans.symbols <- mapIds(org.Hs.eg.db, keys = unique(trans_eQTL$ENSG),
                            column = "SYMBOL", keytype = "ENSEMBL",
                            multiVals = "first") %>% na.omit

trans_eQTL <- trans_eQTL %>%
    mutate(GeneSymbol = trans.symbols[ENSG])

message("[4/7] Selecting expressed genes trans-eQTLs and calculating beta")
setkey(trans_eQTL, GeneSymbol)
trans_eQTL <- trans_eQTL[SNP %in% cis_eSNPs$SNP & GeneSymbol %in% perturb_summary_stats$effect]
message("  ✔ Matching trans-eQTLs: ", nrow(trans_eQTL))

message("[5/7] Merging cis and trans-eQTLs of the lead cis-eSNPs of perturbed genes")
merged.QTL <- merge(
  cis_eSNPs,
  trans_eQTL,
  allow.cartesian = TRUE,
  by = "SNP", suffixes = c(".perturb", ".effect")
) %>% dplyr::rename(
  perturb = GeneSymbol.perturb,
  effect = GeneSymbol.effect,
)

message("  ✔ eQTL pairs: ", nrow(merged.QTL))

message("[6/7] Selecting lead cis-eSNP of perturbed gene")

message("[7/7] Merging perturbation and eQTL pairs")
dat <- merge(perturb_summary_stats, merged.QTL, by = c("perturb", "effect"))
message("  ✔ Final perturbation/eQTL data rows: ", nrow(dat))

message("✔ Writing final outputs...")
fwrite(dat, file = output_file)
message("🎉 Done.")
