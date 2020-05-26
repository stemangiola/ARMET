#!/bin/sh
cd $PBS_O_WORKDIR;
# rm -rf /tmp/*
module load R/3.6.1;
Rscript dev/run_armet_on_cancer_type.R $1
