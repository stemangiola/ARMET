#!/bin/sh
#cd $PBS_O_WORKDIR;

module load R/4.0.2;
Rscript dev/run_armet_on_cancer_type.R $1
