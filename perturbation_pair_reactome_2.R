# Filter perturbation-effect pairs to only pairs that participate in known pathways in Reactome
# Load necessary libraries
library(ReactomePA)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(reactome.db)
library(tidyverse)
library(here)

# Read the gene interaction pairs
gene_interactions <- read_csv(here("data/perturb/pertubration_pairs.csv"))

# Check structure
print(head(gene_interactions))

# Convert gene symbols to UniProt IDs
all_genes <- unique(c(gene_interactions$perturbation, gene_interactions$effect))

gene_uniprot <- mapIds(org.Hs.eg.db, keys = all_genes, 
                      column = "UNIPROT", keytype = "SYMBOL", 
                      multiVals = "first")

# Remove NAs
valid_uniprot <- na.omit(gene_uniprot)

# Read UniProt to Reactome Pathway ID mapping (the lowest level pathway diagram or subset of the pathway)
# From: https://reactome.org/download-data ; https://reactome.org/download/current/UniProt2Reactome.txt 
uniprot2reactome <- read_tsv(here("data/Reactome/UniProt2Reactome.txt"),
			     col_names=c("uniprot_id",
					 "reactome_pathway_id",
					 "reactome_url",
					 "reactome_pathway_name",
					 "evidence_code",
					 "species")) %>% filter(species == "Homo sapiens")

# Filter to only those proteins that are in the perturbation pairs
uniprot2reactome <- uniprot2reactome %>% filter(uniprot_id %in% valid_uniprot) 

# group by uniprot ID to get lists of pathways per gene
uniprot2reactome <- uniprot2reactome %>% group_by(uniprot_id) %>%
	summarize(reactome_pathway_ids = list(reactome_pathway_id),
		  reactome_pathway_names = list(reactome_pathway_name),
		  .groups = 'drop'  # To drop the grouping structure after summarizing
		  )

gene_interactions$perturbation_uniprot <- gene_uniprot[gene_interactions$perturbation]
gene_interactions$effect_uniprot <- gene_uniprot[gene_interactions$effect]

# Drop rows where the UniProt ID is NA
gene_interactions <- gene_interactions %>% na.omit

# Merge gene_interactions with uniprot2reactome to get pathway names for perturbation and effect
gene_interactions_with_pathways <- gene_interactions %>%
  left_join(uniprot2reactome, by = c("perturbation_uniprot" = "uniprot_id")) %>%
  left_join(uniprot2reactome, by = c("effect_uniprot" = "uniprot_id"), suffix = c("_perturbation", "_effect"))

# Function to find the intersection of pathway names
find_intersection <- function(pathways_perturbation, pathways_effect) {
  intersect(pathways_perturbation, pathways_effect)
}

# Add intersecting pathways and filter out rows with no intersection
gene_interactions_filtered <- gene_interactions_with_pathways %>%
  mutate(
    intersecting_pathways = mapply(find_intersection, reactome_pathway_names_perturbation, reactome_pathway_names_effect)
  ) %>%
  filter(lengths(intersecting_pathways) > 0)  # Keep only rows where there is an intersection

gene_interactions_pathways <- gene_interactions_filtered %>%
  select(effect, perturbation, intersecting_pathways) %>%
  unnest(intersecting_pathways)


###### GET THE EFFECT SIZES OF THOSE PATHWAYS

perturb_DEG_dir <- here("data/perturb/DEGs/Jurkat")

gene_interactions_pathways <- gene_interactions_pathways %>%
  mutate(DEG_file_path = file.path(perturb_DEG_dir, paste0(perturbation, ".tsv")),
         effect_ensg = mapIds(org.Hs.eg.db, keys = effect,
                              column = "ENSEMBL", keytype = "SYMBOL",
                              multiVals = "first"))

# Read all DEGs files into one large data frame
DEGs_list <- lapply(gene_interactions_pathways$DEG_file_path %>% unique, function(path) {
  df <- read.table(path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  df$effect_ensg <- rownames(df)  # Add row names as a column for merging

  # Extract filename without extension
  perturbation_value <- sub("\\.tsv$", "", basename(path))
  df$perturbation <- perturbation_value  # Assign it to the entire df

  return(df)
})

# Combine all DEGs into a single data frame
DEGs_combined <- do.call(rbind, DEGs_list)

# Merge with gene_interactions_pathways
gene_interactions_pathways <- merge(
  gene_interactions_pathways, 
  DEGs_combined[, c("effect_ensg", "p_val_adj", "avg_log2FC", "perturbation")], 
  by = c("effect_ensg", "perturbation"), 
  all.x = TRUE
)

### GTEx coexpression

gtex_coexp <- read_tsv(here("data/GTEx_coexpression/WholeBlood_TWNs.out.txt")) %>% filter(`Edge type` == 1)
