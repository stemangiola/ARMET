
test_that("check data set",{
	
	# Test fo anylvel that all gnes are in all cel tpes
	expect_lte(
		ARMET::ARMET_ref %>%
			distinct(level, symbol, `Cell type category`) %>%
			count(level, symbol) %>%
			distinct(level, n) %>%
			nrow,
		4
	)
	
	# Test dataset house keeping should be in all levels
	expect_equal(
		ARMET::ARMET_ref %>% distinct(`house keeping`, level) %>% nrow,
		8
	)
	
	# Test dataset house keeping should be in all samples
	expect_equal(
		ARMET::ARMET_ref %>% distinct(`house keeping`, sample) %>% count(`house keeping`) %>% distinct(n) %>% nrow,
		1
	)
	
	# there should be howse keeping
	expect_equal(
		ARMET::ARMET_ref %>% distinct(`house keeping`) %>% nrow,
		2
	)
	
	# there should be howse keeping
	expect_equal(
		ARMET::ARMET_ref %>% filter(level==3) %>% pull(`Cell type category`) %>% unique %>% intersect(c("endothelial", "epithelial", "fibroblast")) %>% length,
		3
	)
	
	# test if any count is NA
	expect_equal(
		ARMET::ARMET_ref %>% filter(count %>% is.na) %>% nrow,
		0
	)
	
	
	
})

library(tidyverse)
library(ARMET)
my_mix =
	ARMET_ref %>% 
	inner_join( (.) %>% distinct(sample) %>% slice(1:2)) %>% select(-level) %>%
	
	# complete
	mutate_if(is.factor, as.character) %>% 
	tidyr::complete(sample, symbol, fill = list(count = 0)) %>%
	mutate(count = as.integer(count))

test_that("check simple run NO hierarchy",{
	
	armet_obj_NO_hierarchy =
		
		# Read data
		my_mix |>
		nest(data = -sample) %>%
		mutate(factor_of_interest = c(0,1)) %>%
		unnest(data) %>%
		
		# Format
		setup_convolved_lm_NON_hierarchical(
			~ factor_of_interest,
			.sample = sample,
			.transcript = symbol,
			.abundance = count,
			reference = 
				readRDS("/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/mangiola.s/ARMET_dev/dev/TCGA_makeflow_pipeline/ref_jian_3_optimisations.rds") %>%
				filter(level==1)
		)
	
	armet_estimate =
		armet_obj_NO_hierarchy |>
		estimate_convoluted_lm() 
	
	
})


test_that("check simple run",{
	
	armet_obj_hierarchical =
		
		# Read data
		my_mix |>
		nest(data = -sample) %>%
		mutate(factor_of_interest = c(0,1)) %>%
		unnest(data) %>%
		
		# Format
		ARMET:::setup_convolved_lm_hierarchical(
			~ factor_of_interest,
			.sample = sample,
			.transcript = symbol,
			.abundance = count,
			reference = readRDS("/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/mangiola.s/ARMET_dev/dev/TCGA_makeflow_pipeline/ref_jian_3_optimisations.rds")
		)
	
	armet_estimate =
		armet_obj_hierarchical |>
		estimate_convoluted_lm_1() 
	
})

test_that("check nk dataset run",{
	
	`%$%` = magrittr::`%$%`
	
	
	armet_obj =
		
		# Read data
		ARMET_ref %>%
		inner_join( (.) %>% dplyr::filter(`Cell type category` == "nk_primed") %>% distinct(sample) %>% slice(0:1)) %>%
		dplyr::select(-level) %>%
		
		mutate(count = as.numeric(count)) %>%
		
		# complete
		mutate_if(is.factor, as.character) %>% 
		tidyr::complete(sample, symbol, fill = list(count = 0)) %>%
		mutate(count = as.integer(count)) |>
		
		# Format
		setup_convolved_lm_hierarchical(
			~ 1,
			.sample = sample,
			.transcript = symbol,
			.abundance = count,
			reference = readRDS("/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/mangiola.s/ARMET_dev/dev/TCGA_makeflow_pipeline/ref_jian_3_optimisations.rds")
		) |>
		estimate_convoluted_lm_1() %>%
		estimate_convoluted_lm_2() %>%
		estimate_convoluted_lm_3()  %>%
		dplyr::filter(level==3) %>%
		dplyr::filter(`Cell type category` == "nk_primed") %>%
		select(-.variable) %>%
		pull(.value) %>% 
		mean %>%
		expect_gt(0.95)
	
})


# test_that("Check accuracy N52",{
# 
# 	N52_ARMET_T =
# 		filter(	readRDS("dev/N52.rds"),	ct == "T") %>%
# 		mutate(count = as.integer(count)) %>%
# 		ARMET_tc(
# 			.sample = sample,
# 			.transcript = symbol,
# 			.abundance = count
# 		) %>%
# 		ARMET_tc_continue(2) %>%
# 		ARMET_tc_continue(3)
# 
# 	expect_gt(
# 		N52_ARMET_T$proportions %>% filter(`Cell type category` == "immune_cell") %>% select(-.variable) %>%  summarise(.value %>% min),
# 		0.93
# 	)
# 
# 	expect_gt(
# 		N52_ARMET_T$proportions %>% filter(`Cell type category` == "t_cell") %>% select(-.variable) %>% summarise(.value %>% mean),
# 		0.79
# 	)
# 
# 	expect_lt(
# 		N52_ARMET_T$proportions %>% select(-.variable) %>% filter(!converged) %>% nrow,
# 		3
# 	)
# 
# })


# test_that("Simulated data and plot",{
# 
# res =
# 	readRDS("dev/noiseless_mix.rds") %>%
# 	gather(transcript, count, -sample, -covariate_2) %>%
# 	ARMET_tc(
# 		~ covariate_2,
# 		sample,
# 		transcript,
# 		count,
# 		iterations = 600,
# 		sampling_iterations = 400
# 	)  %>%
# 	ARMET_tc_continue(2) %>%
# 	ARMET_tc_continue(3)
# 
# expect_equal(
# 	c("b_cell" ,  "mast_cell" , "b_memory",   "neutrophil") %in%
# 	(res %>%
# 		test_differential_composition() %>%
# 		filter(significant) %>%
# 		pull(`Cell type category`) %>%
# 		unique
# 	) %>%
# 		all(),
# 	TRUE
# )
# 
# 
# # Test plotting
# res %>%
# 	test_differential_composition() %>%
# 	plot_polar()
# 
# })

# test_that("censoring",{
# 
# 	mix = readRDS("dev/mix_for_package_test.rds")
# 
# 	res =
# 		mix %>%
# 		mutate(`count mix` = as.integer(`count mix`)) %>%
# 		mutate(run = as.character(run)) %>%
# 		select(-level) %>%
# 
# 		ARMET_tc(
# 			 ~ censored(days, alive),
# 			run,
# 			symbol,
# 			`count mix`,
# 			iterations = 600,
# 			sampling_iterations = 400,
# 			prior_survival_time = mix %>% distinct(run, real_days) %>% pull(real_days)
# 
# 		) %>%
# 		ARMET_tc_continue(2) %>%
# 		ARMET_tc_continue(3)
# 
# 	expect_equal(
# 		c("immune_cell", "b_cell" ,     "b_memory" ,   "b_naive" ) %in%
# 			(res %>%
# 			 	test_differential_composition() %>%
# 			 	filter(significant) %>%
# 			 	pull(`Cell type category`) %>%
# 			 	unique
# 			) %>%
# 			all(),
# 		TRUE
# 	)
# 
# })



