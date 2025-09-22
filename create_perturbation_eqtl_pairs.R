library(tidyverse)
library(glue)
library(here)
library(org.Hs.eg.db)
library(data.table)
library(argparse)
source(here("src/functions_create_perturb_df.R"))
source(here("src/functions_ashr.R"))

setDTthreads(16)

parser <- ArgumentParser()
parser$add_argument("--cells", type = "character", help="Cell type: K562_essential, K562, Jurkat, RPE1", required=T)
parser$add_argument("--ash", action = "store_true", help="Run adaptive shrinkage (ashr)")
parser$add_argument("--bulk", action = "store_true", help="Pseudobulk?")
parser$add_argument("--efficient", action = "store_true", help="Filter efficient perturbations only (>70% perturb)")
parser$add_argument("--unfiltered", action = "store_true", help="Don't filter perturbations-effect pairs for significance")
args <- parser$parse_args()

cells <- args$cells
ash <- args$ash
bulk <- args$bulk
unfiltered <- args$unfiltered
efficient <- args$efficient

ash_prefix = if(ash) {"ash_"} else {""}
efficient_suffix = if(efficient) {"_eff"} else {""}
unfiltered_suffix = if(unfiltered) {"_all"} else {""}
###### PROCESS PERTURB-SEQ DEGS

bulk_prefix <- if(bulk) "bulk_" else ""
perturb_DEG_dir <- here(glue("data/perturb/DEGs/{bulk_prefix}{cells}"))

DEG_files <- list.files(perturb_DEG_dir, full.names = TRUE)

perturbation_pairs_file <- here(glue("data/perturb/pairs/{bulk_prefix}{cells}/{ash_prefix}{bulk_prefix}perturbation_pairs{unfiltered_suffix}{efficient_suffix}.csv"))
QC <- fread(here(glue("data/perturb/QC/{bulk_prefix}{cells}/gene_QC.csv"))) %>%
  filter(n_expr_trt > 7, n_expr_ctrl > 7)


if(!file.exists(perturbation_pairs_file)) {
  if(ash) {
    perturbation_effect_df <- apply_ash(DEG_files, cores=16, significant=!unfiltered, bulk=bulk)
  } else {
    results <- extract_DEGs(DEG_files,
                            cores=16,
                            significant=!unfiltered,
                            bulk=bulk,
                            efficient = efficient)
  }
  
  perturbation_effect_df <- results$perturbation_effect_df
  filter_summary <- results$filter_summary
  
  # keep <- paste(perturbation_effect_df$perturbation, perturbation_effect_df$effect) %in%
  #   paste(QC$perturbation, QC$effect)
  # perturbation_effect_df <- perturbation_effect_df[keep, ]
  # 
  # filter_summary$QC_pairs <- nrow(perturbation_effect_df)
  
  
  # Write valid perturbation-effect pairs
  fwrite(perturbation_effect_df, here(glue("data/perturb/pairs/{bulk_prefix}{cells}/{ash_prefix}{bulk_prefix}perturbation_pairs{unfiltered_suffix}{efficient_suffix}.csv")))
  
  # Write DEGs filter summary
  
  saveRDS(filter_summary, here(glue("data/perturb/pairs/{bulk_prefix}{cells}/filter_summary_perturbation_pairs.RDS"))) 
} else {
  print("Skipping perturbation pairs. Already exist")
  filter_summary <- readRDS(here(glue("data/perturb/pairs/{bulk_prefix}{cells}/filter_summary_perturbation_pairs.RDS")))
  perturbation_effect_df <- fread(perturbation_pairs_file)
}

####### PROCESS EQTLS
### -- Load and filter cis-eQTLs -- ###
message("[1/9] Reading cis-eQTLs")
cis_eQTL <- fread(here("data/summary_stats/cis_eQTLs_eQTLgen.txt"),
                  sep = "\t",
                  select = c("SNP", "Pvalue", "FDR", "GeneSymbol", "Zscore", "NrSamples", "AssessedAllele", "OtherAllele"))
message("  ✔ cis-eQTLs loaded: ", nrow(cis_eQTL))

filter_summary$total_cis_eQTL_SNPs <- uniqueN(cis_eQTL$SNP)
filter_summary$total_cis_eQTL_genes <- uniqueN(cis_eQTL$GeneSymbol)

message("[2/9] Filtering for perturbed genes")
setkey(cis_eQTL, GeneSymbol)
cis_eQTL.perturbations.all <- cis_eQTL[GeneSymbol %in% perturbation_effect_df$perturbation]
message("  ✔ Perturbation cis-eQTLs found: ", nrow(cis_eQTL.perturbations.all))

filter_summary$cis_eQTL_SNPs_perturbed <- uniqueN(cis_eQTL.perturbations.all$SNP)
filter_summary$cis_eQTL_genes_perturbed <- uniqueN(cis_eQTL.perturbations.all$GeneSymbol)

message("[3/9] Filtering significant cis-eQTLs and calculating beta")
cis_eQTL.perturbations.all <- cis_eQTL.perturbations.all[FDR < 0.05]
cis_eQTL.perturbations.all <- calculate_eQTL_beta(cis_eQTL.perturbations.all)

filter_summary$cis_eQTL_SNPs_perturbed_significant <- uniqueN(cis_eQTL.perturbations.all$SNP)
filter_summary$cis_eQTL_genes_perturbed_significant <- uniqueN(cis_eQTL.perturbations.all$GeneSymbol)

message("[4/9] Aggregating cis-eQTLs by SNP")
cis_eQTL.perturbations.aggregated <- cis_eQTL.perturbations.all[
  , lapply(.SD, function(x) list(x)), by = SNP
]

### -- Load and filter trans-eQTLs -- ###
message("[5/9] Reading trans-eQTLs")
trans_eQTL <- fread(here("data/summary_stats/trans_eQTLs_eQTLgen.txt"),
                    sep = "\t",
                    select = c("SNP", "Pvalue", "FDR", "GeneSymbol", "Zscore", "NrSamples", "AssessedAllele", "OtherAllele"))
message("  ✔ trans-eQTLs loaded: ", nrow(trans_eQTL))

filter_summary$total_trans_eQTL_SNPs <- uniqueN(trans_eQTL$SNP)
filter_summary$total_trans_eQTL_genes <- uniqueN(trans_eQTL$GeneSymbol)

message("[6/9] Filtering for differentially expressed effect genes")
setkey(trans_eQTL, GeneSymbol)
trans_eQTL.effect.unfiltered <- trans_eQTL[GeneSymbol %in% perturbation_effect_df$effect]
message("  ✔ Matching trans-eQTLs: ", nrow(trans_eQTL.effect.unfiltered))

filter_summary$trans_eQTL_SNPs_effect <- uniqueN(trans_eQTL.effect.unfiltered$SNP)
filter_summary$trans_eQTL_genes_effect <- uniqueN(trans_eQTL.effect.unfiltered$GeneSymbol)

message("[7/9] Filtering trans-eQTLs with matching cis SNPs and calculating beta")
# Step 1: Keep only trans-eQTLs that share SNPs with significant cis-eQTLs

#trans_eQTL.effect.unfiltered <- trans_eQTL.effect.unfiltered[SNP %in% cis_eQTL.perturbations.all$SNP]

# Step 2: Filter for significant FDR < 0.05 trans-eQTLs

output_path <- here(glue("data/perturb/pairs/{bulk_prefix}{cells}/"))

trans_eQTL.effect.significant <- trans_eQTL.effect.unfiltered[FDR < 0.05]

# Step 3: Calculate beta
trans_eQTL.effect.unfiltered <- calculate_eQTL_beta(trans_eQTL.effect.unfiltered)

filter_summary$total_trans_eQTL_SNPs_effect_significant <- uniqueN(trans_eQTL.effect.significant$SNP) 
filter_summary$total_trans_eQTL_genes_effect_significant <- uniqueN(trans_eQTL.effect.significant$GeneSymbol) 

#filter_summary$trans_eQTL_DEG_matching_cis <- nrow(trans_eQTL.effect.unfiltered)
#filter_summary$trans_eQTL_DEG_matching_cis_genes <- uniqueN(trans_eQTL.effect.unfiltered$GeneSymbol)

#filter_summary$trans_eQTL_DEG_matching_cis_signficiant <- nrow(trans_eQTL.effect.significant)
#filter_summary$trans_eQTL_DEG_matching_cis_significant_genes <- uniqueN(trans_eQTL.effect.significant$GeneSymbol)

message("[8/9] Aggregating trans-eQTLs by SNP")
trans_eQTL.effect.aggregated <- trans_eQTL.effect.unfiltered[
  , lapply(.SD, function(x) list(x)), by = SNP
]

### -- Merge and flatten -- ###
message("[9/9] Merging cis and trans and expanding QTL pairs")
setkey(cis_eQTL.perturbations.aggregated, SNP)
setkey(trans_eQTL.effect.aggregated, SNP)

merged.QTL <- merge(
  cis_eQTL.perturbations.aggregated,
  trans_eQTL.effect.aggregated,
  by = "SNP", suffixes = c(".perturb", ".effect")
)

QTL.pairs <- merged.QTL[, {
  cis <- data.table(
    perturbation = GeneSymbol.perturb[[1]],
    Pvalue.perturb = Pvalue.perturb[[1]],
    FDR.perturb = FDR.perturb[[1]],
    Beta.perturb = Beta.perturb[[1]]
  )

  trans <- data.table(
    effect = GeneSymbol.effect[[1]],
    Pvalue.effect = Pvalue.effect[[1]],
    FDR.effect = FDR.effect[[1]],
    Beta.effect = Beta.effect[[1]]
  )

  # Create Cartesian product of cis × trans
  cj <- CJ(cis_id = 1:nrow(cis), trans_id = 1:nrow(trans))

  cbind(
    cis[cj$cis_id],
    trans[cj$trans_id]
  )
}, by = SNP]


filter_summary$eQTL_pairs <- nrow(QTL.pairs)
filter_summary$eQTL_SNPs <- uniqueN(QTL.pairs$SNP)
filter_summary$eQTL_perturbs <- uniqueN(QTL.pairs$perturbation)
filter_summary$eQTL_effects <- uniqueN(QTL.pairs$effect)

# Save results
output_file <- glue("{ash_prefix}{bulk_prefix}eQTL_pairs{unfiltered_suffix}{efficient_suffix}.csv")

message("✔ Writing final outputs...")
fwrite(QTL.pairs, file = file.path(output_path, output_file))
saveRDS(filter_summary, file = file.path(output_path, "filter_summary.RDS"))
message("🎉 Done.")

