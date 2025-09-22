#!/bin/bash

Rscript create_perturbation_eqtl_pairs.R --cells Jurkat --efficient
Rscript create_perturbation_eqtl_pairs.R --cells HepG2 --efficient
Rscript create_perturbation_eqtl_pairs.R --cells K562_GenomeWide --efficient
Rscript create_perturbation_eqtl_pairs.R --cells RPE1 --efficient
