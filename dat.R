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

track_log <- function(df, label, cell = NA, log_file) {
  n_pairs <- nrow(df)
  n_perturbs <- if ("perturb" %in% names(df)) length(unique(df$perturb)) else NA
  n_effects <- if ("effect" %in% names(df)) length(unique(df$effect)) else NA
  
  entry <- data.frame(
    step = label,
    cell = cell,
    n_pairs = n_pairs,
    n_perturbs = n_perturbs,
    n_effects = n_effects,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  
  # Append or create
  if (!file.exists(log_file)) {
    fwrite(entry, log_file)
  } else {
    fwrite(entry, log_file, append = TRUE)
  }
  
  message(sprintf("[%s | %s] %d pairs | %d perturbs | %d effects", 
                  label, ifelse(is.na(cell), "all", cell), n_pairs, n_perturbs, n_effects))
  
  invisible(df)
}

log_file <- file.path(pairs_dir, "progress_log.csv")
if (file.exists(log_file)) file.remove(log_file)


# Select perturbation-gene pairs that pass QC (> 7 expressing cells in perturb and control groups)
QC_pairs_pass <- mclapply(names(gene_QC), function(cell) {
  QC <- gene_QC[[cell]]
  df <- QC %>%
    filter(n_expr_trt > 7,  n_expr_ctrl > 7) %>%
    select(perturbation, gene, effect)
  track_log(df, "QC passed pairs", cell, log_file)
  df
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

      DEGs$base_effect <- base_exp$base_mean_per_cell[match(DEGs$gene, base_exp$gene_id)]
      pert_eff_value <- 1 - 2^(DEGs$log2FoldChange[DEGs$effect == perturbed_gene])
      if (length(pert_eff_value) == 0) pert_eff_value <- 1
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
    }, error = function(e) {
      warning("Error in file ", file, ": ", conditionMessage(e))
      return(NULL)
    })
  }, mc.cores = cores_per_cell, mc.preschedule = FALSE)

  perturb_effects <- Filter(Negate(is.null), perturb_effects)  # remove NULLs
  if (length(perturb_effects) == 0) return(NULL)

  df <- bind_rows(perturb_effects)
  track_log(df, "Processed DEGs (after QC)", cell, log_file)
  df
}, mc.cores = length(cells))

names(results_per_cell) <- cells

perturb_summary_stats <- rbindlist(results_per_cell)

track_log(perturb_summary_stats, "Combined DEGs across all cells", NA, log_file)

####### PROCESS EQTLS
### -- Load and filter cis-eQTLs -- ###
message("[1/7] Reading cis-eQTLs")
cis_eQTL <- fread("/rds/project/rds-csoP2nj6Y6Y/biv22/data/eqtl/cis_eQTLs_eQTLgen.txt",
                  sep = "\t",
                  select = c("SNP", "Pvalue", "FDR", "GeneSymbol", "Zscore", "NrSamples", "AssessedAllele", "OtherAllele"))
message("  ✔ cis-eQTLs loaded: ", nrow(cis_eQTL))

message("[2/7] Selecting significant cis-eSNPs and calculating beta")
setkey(cis_eQTL, GeneSymbol)
cis_eSNPs <- cis_eQTL %>%
  as.data.table %>%
  filter(GeneSymbol %in% perturb_summary_stats$perturb) %>%
  calculate_eQTL_beta %>%
  # group_by(GeneSymbol) %>%
  # filter(Pvalue == min(Pvalue)) %>%
  # filter(abs(Beta) == max(abs(Beta))) %>%
  # ungroup() %>%
  filter(FDR < 0.05)

message("  ✔ Significant perturbation cis-eSNPs found: ", nrow(cis_eSNPs))

message("[3/7] Reading trans-eQTLs")
trans_eQTL <- fread("/rds/project/rds-csoP2nj6Y6Y/biv22/data/eqtl/trans_eQTLs_eQTLgen.txt",
                    sep = "\t",
                    select = c("SNP", "Pvalue", "FDR", "GeneSymbol", "Zscore", "NrSamples", "AssessedAllele", "OtherAllele"))
message("  ✔ trans-eQTLs loaded: ", nrow(trans_eQTL))

message("[4/7] Selecting expressed genes trans-eQTLs and calculating beta")
setkey(trans_eQTL, GeneSymbol)
trans_eQTL <- trans_eQTL[SNP %in% cis_eSNPs$SNP & GeneSymbol %in% perturb_summary_stats$effect] %>% calculate_eQTL_beta
message("  ✔ Matching trans-eQTLs: ", nrow(trans_eQTL))

message("[5/7] Merging cis and trans-eQTLs of the lead cis-eSNPs of perturbed genes")
merged.QTL <- merge(
  cis_eSNPs,
  trans_eQTL,
  allow.cartesian = TRUE,
  by = "SNP", suffixes = c(".perturb", ".effect")
) %>% rename(
  perturb = "GeneSymbol.perturb",
  effect = "GeneSymbol.effect",
) %>% mutate(y = ifelse(FDR.effect < 0.05, 1, 0))

message("  ✔ eQTL pairs: ", nrow(merged.QTL))

message("[6/7] Selecting lead cis-eSNP of perturbed gene")

merged.QTL <- merged.QTL %>%
  group_by(perturb) %>%
  filter(Pvalue.perturb == min(Pvalue.perturb)) %>%
  filter(abs(Beta.perturb) == max(abs(Beta.perturb))) %>%
  ungroup()

message("[7/7] Merging perturbation and eQTL pairs")
dat <- merge(perturb_summary_stats, merged.QTL, by = c("perturb", "effect"))
message("  ✔ Final perturbation/eQTL data rows: ", nrow(dat))

message("✔ Writing final outputs...")
fwrite(dat, file = output_file)
message("🎉 Done.")
