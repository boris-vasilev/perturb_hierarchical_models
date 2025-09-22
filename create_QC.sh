#!/bin/bash
. /etc/profile.d/modules.sh
module load gcc/11

Rscript QC.R --cells Jurkat
Rscript QC.R --cells HepG2
Rscript QC.R --cells K562_GenomeWide
Rscript QC.R --cells RPE1
