context('Test if divergent')

library(tidyverse)
library(ARMET)
# library(data.tree)
# library(foreach)
# library(magrittr)


# Test fo anylvel that all gnes are in all cel tpes
expect_lte(
	ARMET::ARMET_ref %>%
		distinct(level, symbol, `Cell type category`) %>%
		count(level, symbol) %>%
		distinct(level, n) %>%
		nrow,
	3
)

# Test dataset house keeping should be in all levels
expect_equal(
	ARMET::ARMET_ref %>% distinct(`house keeping`, level) %>% nrow,
	6
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


my_mix = ARMET_ref %>% inner_join( (.) %>% distinct(sample) %>% slice(1)) %>% distinct(sample, symbol, count) %>% spread(symbol, count)

# result =
# 	ARMET_tc(
# 		my_mix,
# 		iterations = 50,
# 		sampling_iterations = 5,
# 		n_markers = res$input$n_markers,
# 		full_bayesian  = T,
# 		cores = 2, levels = 2
# 	)

result_fix =
	ARMET_tc(
		my_mix,
		iterations = 50,
		sampling_iterations = 5,
		full_bayesian  = F,
		cores = 2,
		levels = 3
	)

result_fix =
	ARMET_tc(
		my_mix,
		iterations = 50,
		sampling_iterations = 5,
		full_bayesian  = T,
		cores = 2,
		levels = 1
	)

# result_dirichlet =
# 	ARMET_tc(
# 		read_csv("dev/mix_dirichlet.csv") %>%
# 			rename(sample = run) %>%
# 			mutate(`count mix` = `count mix` %>% as.integer) %>%
# 			spread(symbol, `count mix`) %>%
# 			mutate(sample = sprintf("s%s", sample))	,
# 		n_markers = res$input$n_markers,
# 		full_bayesian  = T, do_regression = T,X = X,
# 		cores = 10, levels = 3
# 	)

result_dirichlet_fix =
	ARMET_tc(
		read_csv("/stornext/Home/data/allstaff/m/mangiola.s/PhD/deconvolution/ARMET/dev/mix_dirichlet.csv") %>%
			rename(sample = run) %>%
			mutate(`count mix` = `count mix` %>% as.integer) %>%
			spread(symbol, `count mix`) %>%
			mutate(sample = sprintf("s%s", sample))	,
		n_markers = res$input$n_markers,
		full_bayesian  = F,
		do_regression = T,
		X = X,
		cores = 10,
		levels = 3
	)

# N52

N52_ARMET_T =
	ARMET_tc(
		spread(
			select(
				filter(
					readRDS("/stornext/Home/data/allstaff/m/mangiola.s/PhD/deconvolution/ARMET/dev/N52.rds"),
					 ct == "T"
					),
					sample, symbol, count
			),
		symbol, count
		),
		full_bayesian  = F,
		do_regression = T,
		cores = 20,
		levels = 3
	)

expect_gt(
	N52_ARMET_T$proportions %>% filter(`Cell type category` == "immune_cell") %>% summarise(.value %>% min),
	0.95
)

expect_gt(
	N52_ARMET_T$proportions %>% filter(`Cell type category` == "t_cell") %>% summarise(.value %>% min),
	0.83
)

expect_equal(
	N52_ARMET_T$proportions %>% filter(!converged) %>% nrow,
	0
)

N52_ARMET_E =
	ARMET_tc(
		spread(
			select(
				filter(
					readRDS("/stornext/Home/data/allstaff/m/mangiola.s/PhD/deconvolution/ARMET/dev/N52.rds"),
					ct == "E"
				),
				sample, symbol, count
			),
			symbol, count
		),
		full_bayesian  = F,
		do_regression = T,
		cores = 20,
		levels = 1
	)

expect_gt(
	N52_ARMET_E$proportions %>% filter(`Cell type category` == "epithelial") %>% summarise(.value %>% min),
	0.95
)

expect_equal(
	N52_ARMET_E$proportions %>% filter(!converged) %>% nrow,
	0
)

# Melanoma

mela =
	readRDS("/stornext/Home/data/allstaff/m/mangiola.s/PhD/deconvolution/ARMET/dev/melanoma_2_samples.rds") %>%
	ARMET_tc(
		full_bayesian  = F,
		do_regression = T,
		cores = 30,
		levels = 3
	)

expect_equal(
	mela$proportions %>% filter(!converged) %>% nrow,
	0
)
