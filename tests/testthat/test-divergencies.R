context('Test if divergent')

library(rlang)
library(tidyverse)
library(ARMET)
library(data.tree)
library(foreach)
library(magrittr)

res = readRDS("/stornext/Home/data/allstaff/m/mangiola.s/PhD/deconvolution/ARMET/dev/feature_selection_fifth_iterative_run_fix_skylake_2/markers_17.RData")
load("dev/mix_dirichlet_X.rda")


result =
	ARMET_tc(
		res$input$mix %>% slice(1),
		iterations = 50,
		sampling_iterations = 5,
		n_markers = res$input$n_markers,
		full_bayesian  = T,
		cores = 2, levels = 2
	)

result_fix =
	ARMET_tc(
		res$input$mix %>% slice(1),
		iterations = 50,
		sampling_iterations = 5,
		full_bayesian  = F,
		cores = 2, levels = 3
	)


result_dirichlet =
	ARMET_tc(
		read_csv("dev/mix_dirichlet.csv") %>%
			rename(sample = run) %>%
			mutate(`count mix` = `count mix` %>% as.integer) %>%
			spread(symbol, `count mix`) %>%
			mutate(sample = sprintf("s%s", sample))	,
		n_markers = res$input$n_markers,
		full_bayesian  = T, do_regression = T,X = X,
		cores = 10, levels = 3
	)

result_dirichlet_fix =
	ARMET_tc(
		read_csv("dev/mix_dirichlet.csv") %>%
			rename(sample = run) %>%
			mutate(`count mix` = `count mix` %>% as.integer) %>%
			spread(symbol, `count mix`) %>%
			mutate(sample = sprintf("s%s", sample))	,
		n_markers = res$input$n_markers,
		full_bayesian  = F, do_regression = T,
		X = X,
		cores = 10, levels = 3
	)
