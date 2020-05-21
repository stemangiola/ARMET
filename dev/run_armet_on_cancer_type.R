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

i = "ACC.tcga.harmonized.counts.allgenes.rds"
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
	mutate(alive = PFI.2 == 0) 


res = input_df %>%
#	inner_join((.) %>% distinct(sample) %>% slice(1:5)) %>%      # <----------------------------
ARMET_tc(
	~ censored(PFI.time.2, alive),
	sample,
	transcript, 
	count, 
	levels = 1,
	iterations = 1500,                     # <----------------------------
	sampling_iterations = 500,                # <----------------------------
	prior_survival_time = read_csv("dev/survival_TCGA_curated.csv") %>% filter(PFI.2==1 & !is.na(PFI.time.2) & PFI.time.2 != "#N/A") %>% pull(PFI.time.2) %>% as.numeric
	# , 
	# model = rstan::stan_model("~/PhD/deconvolution/ARMET/inst/stan/ARMET_tc_fix_hierarchical.stan", auto_write = F)
)  %>%
	ARMET_tc_continue(2) %>%
	ARMET_tc_continue(3) %>%
	ARMET_tc_continue(4)


save(res, file=sprintf("dev/armet_%s.rda", i), compress = "gzip")

# res %>% plot_scatter() + geom_text()
# # 


# 

# Density
(res$proportions %>%
		unnest(draws) %>%
		filter(A == 2) %>%
		ggplot(aes(.value, color=`Cell type category`)) +
		geom_density() +
		facet_wrap(~.variable, scale="free_y")
) %>% plotly::ggplotly()
#


