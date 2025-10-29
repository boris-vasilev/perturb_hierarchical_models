library(tidyverse)
library(data.table)
library(parallel)
library(ashr)
library(argparse)
source("functions_create_perturb_df.R")

parser <- ArgumentParser()
parser$add_argument("--cores", type = "numeric", help = "Number of cores", required = TRUE)

args <- parser$parse_args()

n_cores <- args$cores

cells <- c("HepG2", "Jurkat", "K562_essential", "RPE1")

data_dir <- "/rds/project/rds-csoP2nj6Y6Y/biv22/data"
DE_dir <- file.path(data_dir, "perturb")
QC_dir <- file.path(data_dir, "perturb_QC")
pairs_dir <- file.path(data_dir, "pairs")

output_file <- file.path(pairs_dir, "full_dat.csv")

cell_DE_dirs <- setNames(file.path(DE_dir, cells), cells)

cell_DEG_files <- lapply(cell_DE_dirs, function(dir) {list.files(dir, full.names = TRUE)})

gene_QC_files <- setNames(file.path(QC_dir, cells, "gene_QC.csv"), cells)
base_expression_files <- setNames(file.path(QC_dir, cells, "base_expression.csv"), cells)

gene_QC <- mclapply(gene_QC_files, fread, mc.cores = length(cells))
base_expressions <- mclapply(base_expression_files, fread, mc.cores = length(cells))

# Select perturbation-gene pairs that pass QC (> 7 expressing cells in perturb and control groups)
QC_pairs_pass <- mclapply(gene_QC, function(QC) {
  QC %>%
    filter(n_expr_trt > 7,  n_expr_ctrl > 7) %>%
    select(perturbation, gene, effect)
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
    perturbed_gene <- basename(file) %>% sub("\\.tsv$", "", .)

    message("Reading file: ", file)
    # read file
    DEGs <- fread(file)
    required_cols <- c("gene", "log2FoldChange", "lfcSE", "padj")
    if (!all(required_cols %in% colnames(DEGs))) {
      warning("Skipping ", file, " (missing columns)")
      return(NULL)
    }

    # filter by QC pairs
    DEGs <- DEGs %>%
      inner_join(QC_pass %>% filter(perturbation == perturbed_gene),
                  by = "gene")

    if (nrow(DEGs) == 0) return(NULL)

    # shrink estimates
    ash_fit <- tryCatch(
      ashr::ash(beta = DEGs$log2FoldChange, se = DEGs$lfcSE, method = "shrink"),
      error = function(e) NULL
    )

    if (is.null(ash_fit)) return(NULL)
    # add base expression
    DEGs$base_effect <- base_exp$base_mean_per_cell[match(DEGs$gene, base_exp$gene_id)]

    # add perturbation efficiency
    pert_eff_value <- 1 - 2^(DEGs$log2FoldChange[DEGs$effect == perturbed_gene])
    if (length(pert_eff_value) == 0) pert_eff_value <- 1  # 100% efficiency if the perturbed gene is not expressed at all
    DEGs$perturb_eff <- pert_eff_value

    data.frame(
      perturb = perturbed_gene,
      gene = DEGs$gene,
      effect = DEGs$effect,
      base_effect = DEGs$base_effect,
      perturb_eff = DEGs$perturb_eff,
      logFC = DEGs$log2FoldChange,
      padj = DEGs$padj,
      x = ifelse(DEGs$padj < 0.05, 1, 0),
      ash_logFC = ash_fit$result$PosteriorMean,
      ash_svalue = ash_fit$result$svalue,
      ash_x = ifelse(ash_fit$result$svalue < 0.05, 1, 0),
      screen = cell
    )
  }, mc.cores = cores_per_cell)

  # combine results for this screen
  bind_rows(perturb_effects)
}, mc.cores = length(cells))

names(results_per_cell) <- cells

perturb_summary_stats <- rbindlist(results_per_cell)

####### PROCESS EQTLS
### -- Load and filter cis-eQTLs -- ###
message("[1/9] Reading cis-eQTLs")
cis_eQTL <- fread("/rds/project/rds-csoP2nj6Y6Y/biv22/data/eqtl/cis_eQTLs_eQTLgen.txt",
                  sep = "\t",
                  select = c("SNP", "Pvalue", "FDR", "GeneSymbol", "Zscore", "NrSamples", "AssessedAllele", "OtherAllele"))
message("  ✔ cis-eQTLs loaded: ", nrow(cis_eQTL))

message("[2/9] Selecting lead cis-eSNPs and calculating beta")
setkey(cis_eQTL, GeneSymbol)
lead_cis_eSNPs <- cis_eQTL %>%
  as.data.table %>%
  filter(GeneSymbol %in% perturb_summary_stats$perturb) %>%
  calculate_eQTL_beta %>%
  group_by(GeneSymbol) %>%
  filter(Pvalue == min(Pvalue)) %>%
  filter(abs(Beta) == max(abs(Beta))) %>%
  ungroup() %>%
  filter(FDR < 0.05)

message("  ✔ Lead perturbation cis-eSNPs found: ", nrow(lead_cis_eSNPs))

message("[3/9] Reading trans-eQTLs")
trans_eQTL <- fread("/rds/project/rds-csoP2nj6Y6Y/biv22/data/eqtl/trans_eQTLs_eQTLgen.txt",
                    sep = "\t",
                    select = c("SNP", "Pvalue", "FDR", "GeneSymbol", "Zscore", "NrSamples", "AssessedAllele", "OtherAllele"))
message("  ✔ trans-eQTLs loaded: ", nrow(trans_eQTL))

message("[4/9] Selecting expressed genes trans-eQTLs and calculating beta")
setkey(trans_eQTL, GeneSymbol)
trans_eQTL <- trans_eQTL[SNP %in% lead_cis_eSNPs$SNP & GeneSymbol %in% perturb_summary_stats$effect] %>% calculate_eQTL_beta
message("  ✔ Matching trans-eQTLs: ", nrow(trans_eQTL))

message("[5/9] Merging cis and trans-eQTLs of the lead cis-eSNPs of perturbed genes")
merged.QTL <- merge(
  lead_cis_eSNPs,
  trans_eQTL,
  by = "SNP", suffixes = c(".perturb", ".effect")
) %>% rename(
  perturb = "GeneSymbol.perturb",
  effect = "GeneSymbol.effect",
) %>% mutate(y = ifelse(FDR.effect < 0.05, 1, 0))

message("  ✔ eQTL pairs: ", nrow(merged.QTL))

message("[6/6] Merging perturbation and eQTL pairs")
dat <- merge(perturb_summary_stats, merged.QTL, by = c("perturb", "effect"))
message("  ✔ Final perturbation/eQTL data rows: ", nrow(dat))

message("✔ Writing final outputs...")
fwrite(dat, file = output_file)
message("🎉 Done.")
