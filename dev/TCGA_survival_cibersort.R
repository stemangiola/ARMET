#~/unix3XX/third_party_sofware/cctools-7.1.5-x86_64-centos6/bin/makeflow  -T torque -B '-l walltime=148:60:00' --do-not-save-failed-output dev/makefile_run_ARMET_on_TCGA

args = commandArgs(trailingOnly=TRUE)
library(tidyverse)
library(ARMET)
library(furrr)
library(tidybulk)

plan(multicore)

# PFI the best



# read_csv("dev/survival_TCGA_evaluation.csv") %>%
# 	gather(survival_type, evaluation, c("OS" ,  "PFI" ,  "DFI" ,  "DSS")) %>%
# 	filter(evaluation > 0) %>%
# 	count(survival_type)

my_dir = "../"
my_dir = "~/unix3XX/PhD/deconvolution/"   # <----------------------------


# # Build HPC list
# dir(my_dir) %>%
# 	#grep("_Primary_Tumor", ., value = T) %>%
# 	gsub(".csv", "", .) %>%
# 	sort %>%
# 	sprintf("qsub -l nodes=1:ppn=24,mem=30gb,walltime=148:60:00 dev/run_armet_on_cancer_type.sh -F %s", .) %>%
# 	write_lines("dev/HPC_run_ARMET_on_TCGA.sh")

# dir(sprintf("%sTCGA_harmonised",my_dir)) %>%
# 	sort %>%
# 	sprintf("dev/armet_%s.rda:\n\tdev/run_armet_on_cancer_type.sh %s", ., .) %>%
# 	write_lines("dev/makefile_run_ARMET_on_TCGA")

# armet_Adrenocortical_Carcinoma_Primary_Tumor.rda:
# 	dev/run_armet_on_cancer_type.sh Adrenocortical_Carcinoma_Primary_Tumor

i = "MESO.tcga.harmonized.counts.allgenes.rds"
i = args[1]

outliers = c("TCGA-12-3652", "TCGA-02-2485", "TCGA-12-0618", "TCGA-19-1390", "TCGA-15-1444", "TCGA-41-2571", "TCGA-28-2499")

input_df = 
	readRDS(sprintf("%s/TCGA_harmonised/%s", my_dir, i)) %>%
	as_tibble(rownames = "ens") %>%
	gather(sample, count, -ens) %>%
	mutate(count = as.integer(count)) %>%
	tidyr::extract(sample, into = "sample", regex = "([a-zA-Z0-9]+-[a-zA-Z0-9]+-[a-zA-Z0-9]+)") %>%
	
	# Select primary tumour
	inner_join(
		dir(sprintf("%s/TCGA_harmonised_clinical", my_dir), full.names = T) %>% 
			map_dfr(~ .x %>% readRDS %>% distinct(sample, definition))  %>% 
			filter(definition == "Primary solid Tumor") %>%
			tidyr::extract(sample, into = "sample", regex = "([a-zA-Z0-9]+-[a-zA-Z0-9]+-[a-zA-Z0-9]+)") 
	) %>%
	
	ensembl_to_symbol(ens) %>%
	left_join(
		read_csv("dev/survival_TCGA_curated.csv") %>% 
			select(bcr_patient_barcode, type, PFI.2, PFI.time.2) %>%
			mutate(PFI.2 = ifelse(PFI.2 == "#N/A", NA, PFI.2)) %>%
			mutate(PFI.time.2 = ifelse(PFI.time.2 == "#N/A", NA, PFI.time.2)) %>%
			mutate(PFI.2 = as.integer(PFI.2), PFI.time.2 = as.integer(PFI.time.2)), 
		by = c("sample" = "bcr_patient_barcode")
	) %>%
	filter(PFI.time.2 %>% is.na %>% `!`) %>%
	filter(sample %in% outliers %>% `!`) %>%
	#mutate_if(is.character, as.factor) %>%
	
	# Aggregate duplcates
	aggregate_duplicates(sample, transcript, count) %>%
	mutate(alive = PFI.2 == 0) %>%
	
	# Filter 0 time
	filter(PFI.time.2 != 0)

# Correct for 0 prop ##############################
###################################################

#https://www.rdocumentation.org/packages/DirichletReg/versions/0.3-0/topics/DR_data
compless_0_proportions = function(proportion){
	proportion = matrix(proportion, ncol=1)
	
	( proportion*(nrow(proportion)-1) + (1/ncol(proportion)) ) / nrow(proportion)
}



res = input_df %>%
	tidybulk::deconvolve_cellularity(sample, transcript, count) 


x =	res %>%
		nanny::subset(sample) %>%
		pivot_longer(names_prefix = "cibersort: ", cols = starts_with("cibersort"), names_to = "cell_type", values_to = "proportion") %>%
		nest(data = -cell_type) %>%
		
		# Test survival
		mutate(surv_test = map(data, ~
													 
													 	.x %>% 
				when(filter(., proportion!=0) %>% nrow %>% `>` (0) ~ 
				(.) %>% 
					mutate(dead = if_else(alive, 0, 1)) %>%
					mutate(min_proportion = min(proportion[proportion!=0])) %>%
					mutate(proportion = if_else(proportion==0, min_proportion, proportion)) %>%
					mutate(proportion = proportion  %>% boot::logit()) %>%
					coxph(Surv(PFI.time.2, dead) ~ proportion, .) %>%
					broom::tidy()
												)	 
													 
		)) %>%
		unnest(surv_test)

