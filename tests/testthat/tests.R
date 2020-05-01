context('Test if divergent')

library(tidyverse)
library(ARMET)
# library(data.tree)
# library(foreach)
# library(magrittr)

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

my_mix = ARMET_ref %>% inner_join( (.) %>% distinct(sample) %>% slice(1))


test_that("check simple run",{

result_fix =
	ARMET_tc(
		my_mix,
		.sample = sample,
		.transcript = symbol,
		.abundance = count,
		iterations = 50,
		sampling_iterations = 5,
		cores = 2
	)


})

test_that("check nk dataset run",{

	`%$%` = magrittr::`%$%`

	result_nk_fix =
		ARMET_ref %>%
		inner_join( (.) %>% filter(`Cell type category` == "nk_primed") %>% distinct(sample) %>% slice(1)) %>%
		ARMET_tc(
			.sample = sample,
			.transcript = symbol,
			.abundance = count,
			cores = 2
		)  %>%
		ARMET_tc_continue(2) %>%
		ARMET_tc_continue(3) %$%
		proportions %>%
		filter(level==3) %>%
		filter(`Cell type category` == "nk_primed") 
	
	expect_gt(
		result_nk_fix %>%
			select(-.variable) %>%
			unnest(proportions) %>%
			pull(.value),
		0.95
	)
	
	

	})


test_that("Check accuracy N52",{

	N52_ARMET_T =
		filter(
			readRDS("/stornext/Home/data/allstaff/m/mangiola.s/PhD/deconvolution/ARMET/dev/N52.rds"),
			ct == "T"
		) %>%
		ARMET_tc(
			.sample = sample,
			.transcript = symbol,
			.abundance = count,
			cores = 20,
			do_regression = T
		) %>%
		ARMET_tc_continue(2) %>%
		ARMET_tc_continue(3)

	expect_gt(
		N52_ARMET_T$proportions %>% filter(`Cell type category` == "immune_cell") %>% summarise(.value %>% min),
		0.95
	)

	expect_gt(
		N52_ARMET_T$proportions %>% filter(`Cell type category` == "t_cell") %>% summarise(.value %>% min),
		0.7
	)

	expect_lt(
		N52_ARMET_T$proportions %>% filter(!converged) %>% nrow,
		7
	)

})




test_that("Simulated data and plot",{

res =
	readRDS("dev/noiseless_mix.rds") %>%
	gather(transcript, count, -sample, -covariate_2) %>%
	ARMET_tc(
		~ covariate_2,
		sample,
		transcript,
		count,
		do_regression = T,
		iterations = 600,
		sampling_iterations = 400
	)  %>%
	ARMET_tc_continue(2) %>%
	ARMET_tc_continue(3)

expect_equal(
	c("b_cell" ,  "mast_cell" , "b_memory",   "neutrophil") %in%
	(res %>%
		test_differential_composition() %>%
		filter(significant) %>%
		pull(`Cell type category`) %>%
		unique
	) %>%
		all(),
	TRUE
)


# Test plotting
res %>%
	test_differential_composition() %>%
	plot_polar()

})





test_that("censoring",{

	res =
		readRDS("dev/noiseless_mix.rds") %>%
		mutate(alive = as.integer(sample(c(0,1), n(), replace=T))) %>%
		mutate(covariate_2 = covariate_2 + 1) %>%
		gather(transcript, count, -sample, -covariate_2, -alive) %>%
		ARMET_tc(
			 ~ censored(covariate_2, alive),
			sample,
			transcript,
			count,
			do_regression = T,
			iterations = 600,
			sampling_iterations = 400,
			prior_survival_time = sample(1, 2, 30, replace = T)
		
		) %>%
		ARMET_tc_continue(2) %>%
		ARMET_tc_continue(3)

	expect_equal(
		c("b_cell" ,  "mast_cell" , "b_memory",   "granulocyte") %in%
			(res %>%
			 	test_differential_composition() %>%
			 	filter(significant) %>%
			 	pull(`Cell type category`) %>%
			 	unique
			) %>%
			all(),
		TRUE
	)

})



