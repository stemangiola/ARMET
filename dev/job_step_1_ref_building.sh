cd $PBS_O_WORKDIR
module load R/3.6.0
Rscript dev/step_1_ref_building.R $1 > log_$1.txt
