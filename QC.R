library(Seurat)
library(sceasy)
library(reticulate)
library(tidyverse)
library(argparse)
source("functions_ashr.R")
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflicts_prefer(base::setdiff)

use_condaenv('env_nf')

parser <- ArgumentParser()
parser$add_argument("--cells", type = "character", help = "Cell type/line", required = TRUE)
args <- parser$parse_args()
cells <- args$cells

read_perturb_h5ad <- function(cells) {
  if(cells == "Jurkat") {
    h5ad_file_name <- "jurkat_raw_singlecell_01.h5ad"
  } else if(cells == "K562_GenomeWide") {
    h5ad_file_name <- "K562_gwps_raw_singlecell.h5ad"
  } else if(cells == "HepG2") {
    h5ad_file_name <- "hepg2_raw_singlecell_01.h5ad"
  } else if(cells == "RPE1") {
    h5ad_file_name <- "rpe1_raw_singlecell_01.h5ad"
  }
  return(
    sceasy::convertFormat(file.path("/rds/project/rds-csoP2nj6Y6Y/biv22/perturb_pseudobulk_de/data", h5ad_file_name), from="anndata", to="seurat")
  )
}

seurat_obj <- read_perturb_h5ad(cells)

DEGs <- list.files(file.path("/rds/project/rds-csoP2nj6Y6Y/biv22/data/perturb", cells), full.names = T)

pairs <- extract_DEGs(DEGs, cores = 32, bulk=T)$perturbation_effect_df %>%
  select(perturbation, effect_ensg)

gene_QC_metrics <- function(seurat_object, pairs, expression_threshold = 0, min_cells = 7) {
  expr_mat <- GetAssayData(seurat_object, slot = "data")
  meta <- seurat_object@meta.data
  
  # Filter pairs to valid ones
  valid_genes <- rownames(expr_mat)
  valid_perts <- unique(meta$gene)
  pairs <- pairs %>%
    filter(effect_ensg %in% valid_genes, perturbation %in% valid_perts)
  
  cat("Total valid pairs:", nrow(pairs), "\n")
  
  # Get control expression
  ctrl_cells <- rownames(meta)[meta$gene == "non-targeting"]
  ctrl_expr_mat <- expr_mat[, ctrl_cells, drop = FALSE]
  ctrl_counts <- Matrix::rowSums(ctrl_expr_mat > expression_threshold)
  ctrl_df <- tibble(gene = names(ctrl_counts), n_expr_ctrl = ctrl_counts)
  
  # Group treatment cells by perturbation
  pert_cell_index <- split(rownames(meta), meta$gene)
  pert_cell_index <- pert_cell_index[names(pert_cell_index) != "non-targeting"]
  
  trt_names <- names(pert_cell_index)
  total_perts <- length(trt_names)
  
  # Process each perturbation
  trt_counts_list <- mclapply(seq_along(trt_names), function(i) {
    pert <- trt_names[i]
    cat(sprintf("Processing %d / %d: %s\n", i, total_perts, pert))
    
    pert_cells <- pert_cell_index[[pert]]
    trt_expr_mat <- expr_mat[, pert_cells, drop = FALSE]
    counts <- Matrix::rowSums(trt_expr_mat > expression_threshold)
    
    tibble(gene = names(counts), perturbation = pert, n_expr_trt = counts)
  }, mc.cores = 32, mc.preschedule = FALSE)
  
  trt_counts_df <- bind_rows(trt_counts_list)
  
  # Join with pair list and control counts
  df <- pairs %>%
    rename(gene = effect_ensg) %>%
    left_join(trt_counts_df, by = c("perturbation", "gene")) %>%
    left_join(ctrl_df, by = "gene") %>%
    filter(!is.na(n_expr_trt), !is.na(n_expr_ctrl)) %>%
    mutate(QC_pass = ifelse(n_expr_trt >= min_cells & n_expr_ctrl >= min_cells, "Pass", "Fail"))
  
  return(df)
}

gene_QC <- gene_QC_metrics(seurat_obj, pairs)


DRGs <- read_csv(file.path("/rds/project/rds-csoP2nj6Y6Y/biv22/data/perturb/pairs", cells, "perturbation_pairs_eff.csv"))

map_ensg_to_symbol <- function(gene_QC) {
  gene.ensg <- gene_QC$gene %>% unique
  
  gene.symbols <- mapIds(org.Hs.eg.db, keys = gene.ensg,
                            column = "SYMBOL", keytype = "ENSEMBL",
                            multiVals = "first")
  
  gene.symbols <- gene.symbols %>% na.omit
  
  gene_QC <- gene_QC %>%
    mutate(effect = gene.symbols[gene])
  
  return(gene_QC)
}

gene_QC <- gene_QC %>% map_ensg_to_symbol %>% filter(perturbation != effect)


gene_QC <- left_join(
  gene_QC,
  DRGs,
  by = c("perturbation", "effect")
) %>%
  select(perturbation, effect, n_expr_trt, n_expr_ctrl, avg_log2FC, p_val_adj)

write_csv(gene_QC, file.path("/rds/project/rds-csoP2nj6Y6Y/biv22/data/perturb_QC", cells, "gene_QC.csv"))
