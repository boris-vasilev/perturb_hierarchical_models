#!/bin/bash

Rscript continuous_dat.R --cells Jurkat --efficient --response ratio
Rscript continuous_dat.R --cells HepG2 --efficient --response ratio
Rscript continuous_dat.R --cells K562_GenomeWide --efficient --response ratio
Rscript continuous_dat.R --cells RPE1 --efficient --response ratio

Rscript logistic_dat.R --cells Jurkat --efficient
Rscript logistic_dat.R --cells HepG2 --efficient
Rscript logistic_dat.R --cells K562_GenomeWide --efficient
Rscript logistic_dat.R --cells RPE1 --efficient
