#!/bin/bash

# Rscript create_perturbation_eqtl_pairs.R --cells Jurkat --efficient
# Rscript create_perturbation_eqtl_pairs.R --cells HepG2 --efficient
# Rscript create_perturbation_eqtl_pairs.R --cells K562_GenomeWide --efficient
# Rscript create_perturbation_eqtl_pairs.R --cells K562_essential --efficient
# Rscript create_perturbation_eqtl_pairs.R --cells RPE1 --efficient

# Create unfiltered (any FDR) set of perturbation pairs
Rscript create_perturbation_eqtl_pairs.R --cells Jurkat --efficient --unfiltered
Rscript create_perturbation_eqtl_pairs.R --cells HepG2 --efficient --unfiltered
Rscript create_perturbation_eqtl_pairs.R --cells K562_GenomeWide --efficient --unfiltered
Rscript create_perturbation_eqtl_pairs.R --cells K562_essential --efficient --unfiltered
Rscript create_perturbation_eqtl_pairs.R --cells RPE1 --efficient --unfiltered