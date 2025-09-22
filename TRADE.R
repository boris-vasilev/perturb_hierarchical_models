library(tidyverse)
library(here)
library(argparse)
library(glue)
library(TRADEtools)
library(parallel)


parser <- ArgumentParser()
parser$add_argument("--cells", type = "character", help = "Cell type/line", required = TRUE)
parser$add_argument("--cores", type = "numeric", help = "How many CPUs to run on", default = 2)

args <- parser$parse_args()
cells <- args$cells
cores <- args$cores

run_TRADE_single <- function(file) {
  tryCatch({
  DEGs <- read.table(file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  perturbed_gene <- basename(file) %>% sub("\\.tsv$", "", .)
  TRADE_results <- TRADE(mode = "univariate",
                         results1 = DEGs)
  saveRDS(TRADE_results, file=here(glue("data/perturb/TRADE/{cells}/{perturbed_gene}.rds")))

  return(data.frame(perturbation=perturbed_gene,
		    TI=TRADE_results$distribution_summary$transcriptome_wide_impact,
		    aff_genes=TRADE_results$distribution_summary$Me))
  }, error = function(e) {
    write(glue("TRADE failed for {file}: {e$message}"), file = stderr())
    return(NULL)
  })
}

print(glue("Running TRADE: Cells: {cells}"))

print("Reading perturbation DEGs")
perturb_DEG_dir <- here(glue("data/perturb/DEGs/{cells}"))
DEG_files <- list.files(perturb_DEG_dir, full.names = TRUE)

TRADE_summary <- mclapply(DEG_files, run_TRADE_single, mc.cores = cores)

TRADE_summary <- bind_rows(TRADE_summary)

write_tsv(TRADE_summary, file=here(glue("data/perturb/TRADE/summary_{cells}.tsv")))
