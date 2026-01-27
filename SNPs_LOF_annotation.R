# Get functional consequences of SNPs of eQTLs in perturbation pairs

library(tidyverse)
library(here)
library(biomaRt)
library(data.table)
library(rhdf5)

# Read GWPS pseudobulk expression data
# Not used directly for analysis just to rename the perturbations (row names) of the AD test statistic file from ENSG to symbol
perturb.bulk.k562 <- H5Fopen(here('data/perturb/seurat/K562_gwps_raw_singlecell.h5ad'))

gene_name <- factor(h5read(perturb.bulk.k562, "/var/gene_name"), labels=h5read(perturb.bulk.k562, "/var/__categories/gene_name")) %>% as.character
gene_id <- h5read(perturb.bulk.k562, "/var/gene_id")

expressed_gene_name2symbol_mapping <- data.frame(row.names = gene_id, symbol = gene_name)
expressed_gene_name2symbol_mapping["ENSG00000284024", "symbol"] <- "MSANTD7"

# Read GWPS AD values and fix gene names for row and column names
# There are 11258 perturbations (columns)
# and 5529 sequenced genes (rows)
# The format of the data is a (affected gene X perturbation gene) matrix

perturb_AD <- fread(here('data/summary_stats/andersondarlingpvaluesBHcorrected.csv.gz'))

# Ensure the symbol column is excluded from the matrix portion
gene_names <- expressed_gene_name2symbol_mapping[perturb_AD$V1, "symbol"]
perturb_AD <- perturb_AD %>% dplyr::select(-V1)
colnames(perturb_AD) <- sapply(strsplit(colnames(perturb_AD), "_"), function(x) x[2])
rownames(perturb_AD) <- gene_names

# Create a named list of affected genes for each perturbation gene
perturbations.affected_genes <- lapply(seq_len(ncol(perturb_AD)), function(j) {
  as.character(gene_names[perturb_AD[[j]] < 0.05])  # Use the gene names vector and filter by < 0.05 for the AD statistic
})
names(perturbations.affected_genes) <- colnames(perturb_AD)  # Assign column (perturbation) names to the list names

perturbation_effect_df <- stack(perturbations.affected_genes)
colnames(perturbation_effect_df) <- c("effect", "perturbation")

# Identify perturbations that appear as effects (i.e. they successfully perturb their targets)
valid_perturbations <- perturbation_effect_df %>%
  filter(effect == perturbation) %>%
  pull(perturbation)

# Keep only rows where perturbation is in the valid list
perturbation_effect_df <- perturbation_effect_df %>%
  filter(perturbation %in% valid_perturbations)

# filter out recursive edges
perturbation_effect_df <- perturbation_effect_df %>% filter(perturbation != effect)

####### PROCESS EQTLS
cis_eQTL <- fread(here("data/summary_stats/cis_eQTLs_eQTLgen.tsv.gz"), sep = "\t")

# get the cis-eQTLs of perturbation genes
cis_eQTL.perturbations.all <- cis_eQTL %>% filter(GeneSymbol %in% perturbation_effect_df$perturbation)

# get only the significant cis-eQTLs of perturbations
cis_eQTL.perturbations.significant <- cis_eQTL.perturbations.all %>%
  filter(Pvalue < 5e-8)


###### ALL  SNPs
# VCF header
cat("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", 
    file=here("data/SNP_info/perturbations_eqtl_all.vcf"))

# Prepare the data to write (selecting columns)
vcf_data <- cis_eQTL.perturbations.all[, .(SNPChr, SNPPos, SNP, OtherAllele, AssessedAllele, ".", ".", ".")] %>% distinct

# Write VCF lines using fwrite (faster than write.table)
fwrite(vcf_data, file=here("data/SNP_info/perturbations_eqtl_all.vcf"),
       append=TRUE, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

###### ALL SIGNIFICANT SNPs ON PERTURBATIONS
# VCF header
cat("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", 
    file=here("data/SNP_info/perturbations_eqtl_significant.vcf"))

# Prepare the data to write (selecting columns)
vcf_data <- cis_eQTL.perturbations.significant[, .(SNPChr, SNPPos, SNP, OtherAllele, AssessedAllele, ".", ".", ".")] %>% distinct

# Write VCF lines using fwrite (faster than write.table)
fwrite(vcf_data, file=here("data/SNP_info/perturbations_eqtl_significant.vcf"),
       append=TRUE, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

# VCF file annotated with LOF from snpEff (locally due to mismatch of JDK versions on HPC)
# java -Xmx8g -jar snpEff.jar -v -motif GRCh37.75 vcf/perturbations_eqtl.vcf > vcf/perturbations_eqtl.ann.basic.vcf

# Extracting LOF information from snpEff annotations
library(VariantAnnotation)

annotated_vcf_file <- here("data/SNP_info/perturbations_eqtl_all.ann.basic.vcf")
vcf_data <- readVcf(annotated_vcf_file, "hg19")

# Extract the INFO column
info_data <- info(vcf_data)

# Filter out the variants with no LOF annotations
lof_variants <- vcf_data[!sapply(info_data$LOF, function(x) length(x) == 0), ]

# Extract SNP IDs
snp_ids <- rownames(lof_variants)

# Extract Gene Names from LOF annotations
# The gene name is the first part of the LOF annotation (before the first "|", remove the "(" before the gene name)
genes <- sapply(info(lof_variants)$LOF, function(x) strsplit(as.character(x), "\\|")[[1]][1]) %>% gsub("^\\(", "", .)

# Create a data frame with SNP and associated Gene
snp_gene_pairs <- data.frame(SNP = snp_ids, Gene = genes)

# View the SNP -> Gene pairs with LOF effect
head(snp_gene_pairs)

write_tsv(snp_gene_pairs, here("data/SNP_info/LOF_all_snp_gene_pairs.tsv"))
