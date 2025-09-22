library(tidyverse)
library(here)

# Read perturb-seq pairs
perturbation_pairs <- read_csv(here("data/perturb/pairs/Jurkat/pertubration_pairs.csv"))

# Read eQTL pairs
eQTL_pairs <- read_csv(here("data/perturb/pairs/Jurkat/eQTL_pairs.csv"))


