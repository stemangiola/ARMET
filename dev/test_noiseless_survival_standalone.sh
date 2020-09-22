cd $PBS_O_WORKDIR;
# rm -rf /tmp/*
module load R/3.6.1;
Rscript dev/test_noiseless_survival_standalone.R $1 $2
