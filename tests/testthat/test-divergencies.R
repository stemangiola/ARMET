context('Test if divergent')

library(rlang)
library(tidyverse)
library(ARMET)
library(data.tree)
library(foreach)
library(magrittr)

res = readRDS("/stornext/Home/data/allstaff/m/mangiola.s/PhD/deconvolution/ARMET/dev/feature_selection_fifth_iterative_run_fix_skylake_2/markers_17.RData")
result_full =
	ARMET_tc(
		res$input$mix %>% filter(sample=="10_b_cell natural_killer"),
		iterations = 50,
		sampling_iterations = 5,
		n_markers = res$input$n_markers,
		full_bayesian  = T,
		cores = 10, levels = 2
	)

result_full =
	ARMET_tc(
		res$input$mix %>% filter(sample=="10_b_cell natural_killer"),
		iterations = 50,
		sampling_iterations = 5,
		n_markers = res$input$n_markers,
		full_bayesian  = F,
		cores = 10, levels = 2
	)
