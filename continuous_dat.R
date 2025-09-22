library(tidyverse)
library(here)
library(argparse)
library(glue)
library(data.table)


parser <- ArgumentParser()
parser$add_argument("--cells", type = "character", help = "Cell type/line", required = TRUE)
parser$add_argument("--response", type = "character", help = "Response type: trans-eQTL effect size or eQTL effect ratio", required = TRUE, choices = c("trans", "ratio"))
parser$add_argument("--ash", action = "store_true", help="Run adaptive shrinkage (ashr)")
parser$add_argument("--bulk", action = "store_true", help="Pseudobulk")
parser$add_argument("--efficient", action = "store_true", help="Efficient perturbations >70% only")

args <- parser$parse_args()
cells <- args$cells
response <- args$response


print(glue("Generating data file: Cells: {cells} Response: {response}"))

ash <- if(args$ash) "ash_" else ""
bulk <- if(args$bulk) "bulk_" else ""
eff <- if(args$efficient) "_eff" else ""

# Read perturb-seq pairs
perturbation_pairs <- fread(here(glue("data/perturb/pairs/{bulk}{cells}/{ash}{bulk}perturbation_pairs{eff}.csv")))

# Read eQTL pairs
eQTL_pairs <- fread(here(glue("data/perturb/pairs/{bulk}{cells}/{ash}{bulk}eQTL_pairs{eff}.csv")))

# Identify the lead SNP with the smallest perturbation cis-eQTL p-value
lead_snp_pairs <- eQTL_pairs %>%
  group_by(perturbation) %>%
  filter(Pvalue.perturb == min(Pvalue.perturb)) %>%
  filter(abs(Beta.perturb) == max(abs(Beta.perturb))) %>%
  ungroup()

# Select unique perturb-seq pairs
perturb_effect <- perturbation_pairs %>%
  select(perturbation, effect, avg_log2FC) %>%
  filter(perturbation != effect) %>%
  filter(effect != "") # filter out genes mapped to empty gene symbol

# Select only the pairs with significant trans-eQTLs
cis_trans <- lead_snp_pairs %>%
  filter(FDR.effect < 0.05) %>%
  select(perturbation, effect, Beta.perturb, Beta.effect) %>%
  filter(perturbation != effect)

pairs_df <- cis_trans %>% inner_join(perturb_effect, by=c("perturbation", "effect"))

if(response == "ratio") {
  dat <- with(pairs_df, data.frame(x=avg_log2FC,
                    y=Beta.effect/Beta.perturb,
                    perturb=as.factor(perturbation),
                    effect=as.factor(effect)))
} else {
  dat <- with(pairs_df, data.frame(x=avg_log2FC*Beta.perturb,
                    y=Beta.effect,
                    perturb=as.factor(perturbation),
                    effect=as.factor(effect)))
}

fwrite(dat, here(glue("data/perturb/pairs/{bulk}{cells}/{ash}{bulk}gaussian_{response}_dat{eff}0.05.csv")))
