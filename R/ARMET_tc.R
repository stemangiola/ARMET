
#' ARMET-tc main
#'
#' @description This function calls the stan model.
#'
#'
#' @importFrom tibble tibble
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom dplyr mutate_if
#'
#' @importFrom tidyr spread
#' @importFrom tidyr gather
#' @importFrom tidyr drop_na
#'
#' @importFrom tidybayes gather_samples
#' @importFrom tidybayes median_qi
#'
#' @importFrom foreach foreach
#' @importFrom foreach %do%
#'
#' @import data.tree
#'
#' @param mix A matrix
#' @param my_design A matrix
#' @param cov_to_test A character string
#' @param fully_bayesian A boolean
#' @param is_mix_microarray A boolean
#' @param verbose A boolean
#' @param save_report A boolean
#' @param custom_ref A matrix
#' @param multithread A boolean
#' @param do_debug A boolean
#' @param cell_type_root A character string
#' @param choose_internal_ref A design matrix
#' @param omit_regression A boolean
#' @param save_fit A boolean
#' @param seed An integer
#'
#' @return An ARMET object
#'
#' @export
#'
ARMET_tc = function(
	mix,
	reference = NULL,
	full_bayesian = T,
	verbose =                           F,
	omit_regression =                   F,
	save_fit =                          F,
	seed =                              NULL,
	cores = 14,
	iterations = 300,
	sampling_iterations = 100,
	levels = 1:4,
	full_bayes = T,
	n_markers
){

	input = c(as.list(environment()))


	my_theme =
		theme_bw() +
		theme(
			panel.border = element_blank(),
			axis.line = element_line(),
			panel.grid.major = element_line(size = 0.2),
			panel.grid.minor = element_line(size = 0.1),
			text = element_text(size=12),
			legend.position="bottom",
			aspect.ratio=1,
			axis.text.x = element_text(angle = 90, hjust = 1),
			strip.background = element_blank(),
			axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
			axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
		)

	# Global properties - derived by previous analyses of the whole reference dataset
	sigma_intercept = 1.3420415
	sigma_slope = -0.3386389
	sigma_sigma = 1.1720851
	lambda_mu_mu = 5.612671
	lambda_sigma = 7.131593

	# Set up tree structure

	levels_in_the_tree = 1:4

	tree = 	data.tree::Clone(tree) %>%	{
		# Filter selected levels
		data.tree::Prune(., function(x) x$level <= max(levels_in_the_tree) + 1)
		.
	}

	ct_in_nodes = tree %>% data.tree::ToDataFrameTree("name", "level", "C", "count", "isLeaf") %>% as_tibble %>% arrange(level, C) %>% filter(!isLeaf) %>% pull(count)

	# Get the number of leafs for every level
	ct_in_levels = foreach(l=levels_in_the_tree+1, .combine = c) %do% {

		data.tree::Clone(tree) %>%
			ifelse_pipe(
				(.) %>% data.tree::ToDataFrameTree("level") %>% pull(2) %>% max %>% `>` (l),
				~ {.x;  data.tree::Prune(.x, function(x) x$level <= l);	.x 	}
			)  %>%
			data.tree::Traverse(., filterFun = isLeaf) %>%
			length()
	}

	n_nodes = ct_in_nodes %>% length
	n_levels = ct_in_levels %>% length

	singles_lv2 = tree$Get("C1", filterFun = isLeaf) %>% na.omit %>% as.array
	SLV2 = length(singles_lv2)
	parents_lv2 = tree$Get("C1", filterFun = isNotLeaf) %>% na.omit %>% as.array
	PLV2 = length(parents_lv2)

	singles_lv3 = tree$Get("C2", filterFun = isLeaf) %>% na.omit %>% as.array
	SLV3 = length(singles_lv3)
	parents_lv3 = tree$Get("C2", filterFun = isNotLeaf) %>% na.omit %>% as.array
	PLV3 = length(parents_lv3)

	singles_lv4 = tree$Get("C3", filterFun = isLeaf) %>% na.omit %>% as.array
	SLV4 = length(singles_lv4)
	parents_lv4 = tree$Get("C3", filterFun = isNotLeaf) %>% na.omit %>% as.array
	PLV4 = length(parents_lv4)


	# Print overlap descriptive stats
	#get_overlap_descriptive_stats(mix %>% slice(1) %>% gather(`symbol`, `read count`, -sample), reference)

	# Prepare data frames -
	# For Q query first
	# For G house keeing first
	# For GM level 1 first


	reference_filtered =
		ARMET::ARMET_ref %>%
		left_join(n_markers, by=c("ct1", "ct2")) %>%
		filter_reference(mix) %>%
		select(-ct1, -ct2, -rank, -`n markers`) %>%
		distinct %>%

		# Select cell types in hierarchy
		inner_join(
			tree %>%
				data.tree::ToDataFrameTree("Cell type category", "C") %>%
				as_tibble %>%
				select(-1)

		) %>%

		# # Decrease the number of house keeping used
		# anti_join({
		# 	mdf = (.) %>%
		# 		distinct(symbol, `house keeping`) %>%
		# 		filter(`house keeping`)
		#
		# 	withr::with_seed(	123, 	sample_frac(mdf, 0.8)) %>%
		# 		distinct(symbol)
		# }) %>%

		# Filter on level considered
		filter(level %in% levels)


	df =
		bind_rows(
			# Get reference based on mix genes
			reference_filtered %>% mutate(`query` = FALSE),
			mix %>%
				gather(`symbol`, `read count`, -sample) %>%
				inner_join(reference_filtered %>% distinct(symbol) ) %>%
				left_join( reference_filtered %>% distinct(symbol, `house keeping`) ) %>%
				mutate(`Cell type category` = "query") %>%
				mutate(`query` = TRUE)
		)	%>%

		# Add marker symbol indeces
		left_join(
			(.) %>%
				filter(!`house keeping`) %>%
				distinct(`symbol`) %>%
				mutate(M = 1:n())
		) %>%

		# Add sample indeces
		arrange(!`query`) %>% # query first
		mutate(S = factor(sample, levels = .$sample %>% unique) %>% as.integer) %>%

		# Add query samples indeces
		left_join(
			(.) %>%
				filter(`query`) %>%
				distinct(`sample`) %>%
				mutate(Q = 1:n())
		) %>%

		# Add house keeping into Cell type label
		mutate(`Cell type category` = ifelse(`house keeping`, "house_keeping", `Cell type category`)) %>%
		anti_join(
			(.) %>%
				filter(`house keeping` & !`query`) %>%
				distinct(symbol, level) %>%
				group_by(symbol) %>%
				arrange(level) %>%
				slice(2:max(n(), 2)) %>% # take away house keeping from level 2 above
				ungroup()
		) %>%

		# If house keeping delete level infomation
		mutate(level = ifelse(`house keeping`, NA, level)) %>%

		# Create unique symbol ID
		unite(ct_symbol, c("Cell type category", "symbol"), remove = F) %>%

		# Add gene idx
		left_join(
			(.) %>%
				filter(!`query`) %>%
				distinct(`Cell type category`, ct_symbol, `house keeping`) %>%
				arrange(!`house keeping`, ct_symbol) %>% # House keeping first
				mutate(G = 1:n())
		) %>%
		left_join(
			(.) %>%
				filter(!`house keeping` & !`query`) %>%
				distinct(level, symbol) %>%
				arrange(level, symbol) %>%
				mutate(GM = 1:n()) %>%
				select(-level)
		)

	G = df %>% filter(!`query`) %>% distinct(G) %>% nrow()
	GM = df %>% filter(!`house keeping`) %>% distinct(symbol) %>% nrow()


	# For  reference MPI inference

	shards = cores #* 2

	counts_baseline =
		df %>%

		# Eliminate the query part, not the house keeping of the query
		filter(!`query` | `house keeping`)  %>%

		format_for_MPI(shards)

	S = counts_baseline %>% distinct(sample) %>% nrow()
	N = counts_baseline %>% distinct(idx_MPI, `read count`, `read count MPI row`) %>%  count(idx_MPI) %>% summarise(max(n)) %>% pull(1)
	M = counts_baseline %>% distinct(start, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% max

	y_source =
		df %>%
		filter(`query` & !`house keeping`) %>%
		select(S, Q, `symbol`, `read count`, GM, sample) %>%
		left_join(	df %>% filter(!query) %>% distinct(`symbol`, G, `Cell type category`, level, lambda, sigma_raw, GM, C) ) %>%
		arrange(C, Q, symbol) %>%
		mutate(`Cell type category` = factor(`Cell type category`, unique(`Cell type category`)))

	Q = df %>% filter(`query`) %>% distinct(Q) %>% nrow


	counts_baseline_to_linear =
		counts_baseline %>%
		filter_house_keeping_query_if_fixed(full_bayesian) %>%
		arrange(G, S) %>%
		mutate(counts_idx = 1:n()) %>%
		mutate(S = S %>% as.factor %>% as.integer)

	counts_linear = counts_baseline_to_linear %>%  pull(`read count`)
	G_to_counts_linear = counts_baseline_to_linear %>% pull(G)
	G_linear = G_to_counts_linear
	S_linear = counts_baseline_to_linear %>% pull(S)

	CL = length(counts_linear)
	#G = counts_baseline_to_linear %>%  distinct(G) %>% nrow
	S = counts_baseline_to_linear %>% distinct(S) %>% nrow

	# Counts idx for each level for each level
	counts_idx_lv_NA = counts_baseline_to_linear %>% filter(level %>% is.na) %>% pull(counts_idx)
	CL_NA = counts_idx_lv_NA %>% length
	counts_idx_lv_1 = counts_baseline_to_linear %>% filter(level==1) %>% pull(counts_idx)
	CL_1 = counts_idx_lv_1 %>% length
	counts_idx_lv_2 = counts_baseline_to_linear %>% filter(level==2) %>% pull(counts_idx)
	CL_2 = counts_idx_lv_2 %>% length
	counts_idx_lv_3 = counts_baseline_to_linear %>% filter(level==3) %>% pull(counts_idx)
	CL_3 = counts_idx_lv_3 %>% length
	counts_idx_lv_4 = counts_baseline_to_linear %>% filter(level==4) %>% pull(counts_idx)
	CL_4 = counts_idx_lv_4 %>% length

	# Deconvolution, get G only for markers of each level. Exclude house keeping
	G1_linear = counts_baseline %>% filter(level ==1) %>% distinct(G, GM, C) %>% arrange(GM, C) %>% pull(G)
	G1 = G1_linear %>% length
	G2_linear = counts_baseline %>% filter(level ==2) %>% distinct(G, GM, C) %>% arrange(GM, C) %>% pull(G)
	G2 = G2_linear %>% length
	G3_linear = counts_baseline %>% filter(level ==3) %>% distinct(G, GM, C) %>% arrange(GM, C) %>% pull(G)
	G3 = G3_linear %>% length
	G4_linear = counts_baseline %>% filter(level ==4) %>% distinct(G, GM, C) %>% arrange(GM, C) %>% pull(G)
	G4 = G4_linear %>% length

	# Observed mix counts
	y_linear_1 = y_source %>% filter(level ==1) %>% distinct(GM, Q, S, `read count`) %>% arrange(GM, Q) %>% pull(`read count`)
	y_linear_2 = y_source %>% filter(level ==2) %>% distinct(GM, Q, S, `read count`) %>% arrange(GM, Q) %>% pull(`read count`)
	y_linear_3 = y_source %>% filter(level ==3) %>% distinct(GM, Q, S, `read count`) %>% arrange(GM, Q) %>% pull(`read count`)
	y_linear_4 = y_source %>% filter(level ==4) %>% distinct(GM, Q, S, `read count`) %>% arrange(GM, Q) %>% pull(`read count`)

	# Observed mix samples indexes
	y_linear_S_1 = y_source %>% filter(level ==1) %>% distinct(GM, Q, S, `read count`) %>% arrange(GM, Q) %>% pull(S)
	y_linear_S_2 = y_source %>% filter(level ==2) %>% distinct(GM, Q, S, `read count`) %>% arrange(GM, Q) %>% pull(S)
	y_linear_S_3 = y_source %>% filter(level ==3) %>% distinct(GM, Q, S, `read count`) %>% arrange(GM, Q) %>% pull(S)
	y_linear_S_4 = y_source %>% filter(level ==4) %>% distinct(GM, Q, S, `read count`) %>% arrange(GM, Q) %>% pull(S)

	# Lengths indexes
	Y_1 = y_linear_1 %>% length
	Y_2 = y_linear_2 %>% length
	Y_3 = y_linear_3 %>% length
	Y_4 = y_linear_4 %>% length

	# Non centered
	lambda_mu_prior = c(6.2, 1)
	lambda_sigma_prior =  c( log(3.3) , 1)
	lambda_skew_prior =  c( -2.7, 1)
	sigma_intercept_prior = c( 1.9 , 0.1)

	lambda_log = 	  counts_baseline %>% filter(!query) %>% distinct(G, lambda) %>% arrange(G) %>% pull(lambda)
	sigma_inv_log = counts_baseline %>% filter(!query) %>% distinct(G, sigma_raw) %>% arrange(G) %>% pull(sigma_raw)

	# Linear parallelised
	shards = cores = 8
	shards_in_levels = c(8, 8, 16, 8) %>% `*` (1:4 %in% levels)
	is_level_in = shards_in_levels %>% `>` (0) %>% as.integer
browser()
	weights =
		df %>%
		get_level_lpdf_weights %>%
		arrange(level) %>%
		left_join(tibble(level=1:4, shards = shards_in_levels)) %>%
		uncount(shards) %>%
		pull(weight)


	parse_baseline = function(.data, lv){
		.data %>%
		filter(level==lv) %>%
		distinct(sample, symbol, `Cell type category`, level, `read count`, counts_idx, G, GM, S, `house keeping`) %>%
		left_join( tibble(level=levels, shards = shards_in_levels[levels]) ) %>%
		format_for_MPI_from_linear()
	}

	# Count indexes
	# lv 1

	MPI1 = get_MPI_df(counts_baseline_to_linear, y_source, counts_baseline,1)

	counts_idx_lv_1_MPI = MPI1$counts_idx_lv_MPI
	size_counts_idx_lv_1_MPI = MPI1$size_counts_idx_lv_MPI
	counts_G_lv_1_MPI = MPI1$counts_G_lv_MPI
	size_counts_G_lv_1_MPI = MPI1$size_counts_G_lv_MPI
	counts_G_lv_1_MPI_non_redundant = MPI1$counts_G_lv_MPI_non_redundant
	size_counts_G_lv_1_MPI_non_redundant = MPI1$size_counts_G_lv_MPI_non_redundant
	counts_G_lv_1_MPI_non_redundant_reps = MPI1$counts_G_lv_MPI_non_redundant_reps
	counts_S_lv_1_MPI = MPI1$counts_S_lv_MPI
	size_counts_S_lv_1_MPI = MPI1$size_counts_S_lv_MPI
	y_linear_1_MPI = MPI1$y_linear_MPI
	size_y_linear_1_MPI = MPI1$size_y_linear_MPI
	y_linear_S_1_MPI = MPI1$y_linear_S_MPI
	size_y_linear_S_1_MPI = MPI1$size_y_linear_S_MPI
	G1_linear_MPI = MPI1$G_linear_MPI
	size_G1_linear_MPI = MPI1$size_G_linear_MPI


	MPI2 = get_MPI_df(counts_baseline_to_linear, y_source, counts_baseline,2)

	counts_idx_lv_2_MPI = MPI2$counts_idx_lv_MPI
	size_counts_idx_lv_2_MPI = MPI2$size_counts_idx_lv_MPI
	counts_G_lv_2_MPI = MPI2$counts_G_lv_MPI
	size_counts_G_lv_2_MPI = MPI2$size_counts_G_lv_MPI
	counts_G_lv_2_MPI_non_redundant = MPI2$counts_G_lv_MPI_non_redundant
	size_counts_G_lv_2_MPI_non_redundant = MPI2$size_counts_G_lv_MPI_non_redundant
	counts_G_lv_2_MPI_non_redundant_reps = MPI2$counts_G_lv_MPI_non_redundant_reps
	counts_S_lv_2_MPI = MPI2$counts_S_lv_MPI
	size_counts_S_lv_2_MPI = MPI2$size_counts_S_lv_MPI
	y_linear_2_MPI = MPI2$y_linear_MPI
	size_y_linear_2_MPI = MPI2$size_y_linear_MPI
	y_linear_S_2_MPI = MPI2$y_linear_S_MPI
	size_y_linear_S_2_MPI = MPI2$size_y_linear_S_MPI
	G2_linear_MPI = MPI2$G_linear_MPI
	size_G2_linear_MPI = MPI2$size_G_linear_MPI

	MPI3 = get_MPI_df(counts_baseline_to_linear, y_source, counts_baseline,3)

	counts_idx_lv_3_MPI = MPI3$counts_idx_lv_MPI
	size_counts_idx_lv_3_MPI = MPI3$size_counts_idx_lv_MPI
	counts_G_lv_3_MPI = MPI3$counts_G_lv_MPI
	size_counts_G_lv_3_MPI = MPI3$size_counts_G_lv_MPI
	counts_G_lv_3_MPI_non_redundant = MPI3$counts_G_lv_MPI_non_redundant
	size_counts_G_lv_3_MPI_non_redundant = MPI3$size_counts_G_lv_MPI_non_redundant
	counts_G_lv_3_MPI_non_redundant_reps = MPI3$counts_G_lv_MPI_non_redundant_reps
	counts_S_lv_3_MPI = MPI3$counts_S_lv_MPI
	size_counts_S_lv_3_MPI = MPI3$size_counts_S_lv_MPI
	y_linear_3_MPI = MPI3$y_linear_MPI
	size_y_linear_3_MPI = MPI3$size_y_linear_MPI
	y_linear_S_3_MPI = MPI3$y_linear_S_MPI
	size_y_linear_S_3_MPI = MPI3$size_y_linear_S_MPI
	G3_linear_MPI = MPI3$G_linear_MPI
	size_G3_linear_MPI = MPI3$size_G_linear_MPI

	MPI4 = get_MPI_df(counts_baseline_to_linear, y_source, counts_baseline,4)

	counts_idx_lv_4_MPI = MPI4$counts_idx_lv_MPI
	size_counts_idx_lv_4_MPI = MPI4$size_counts_idx_lv_MPI
	counts_G_lv_4_MPI = MPI4$counts_G_lv_MPI
	size_counts_G_lv_4_MPI = MPI4$size_counts_G_lv_MPI
	counts_G_lv_4_MPI_non_redundant = MPI4$counts_G_lv_MPI_non_redundant
	size_counts_G_lv_4_MPI_non_redundant = MPI4$size_counts_G_lv_MPI_non_redundant
	counts_G_lv_4_MPI_non_redundant_reps = MPI4$counts_G_lv_MPI_non_redundant_reps
	counts_S_lv_4_MPI = MPI4$counts_S_lv_MPI
	size_counts_S_lv_4_MPI = MPI4$size_counts_S_lv_MPI
	y_linear_4_MPI = MPI4$y_linear_MPI
	size_y_linear_4_MPI = MPI4$size_y_linear_MPI
	y_linear_S_4_MPI = MPI4$y_linear_S_MPI
	size_y_linear_S_4_MPI = MPI4$size_y_linear_S_MPI
	G4_linear_MPI = MPI4$G_linear_MPI
	size_G4_linear_MPI = MPI4$size_G_linear_MPI

	# # MODEL
	Sys.setenv("STAN_NUM_THREADS" = cores)

	# library(rstan)
	# fileConn<-file("~/.R/Makevars")
	# writeLines(c( "CXX14FLAGS += -O3","CXX14FLAGS += -DSTAN_THREADS", "CXX14FLAGS += -pthread"), fileConn)
	# close(fileConn)
	# ARMET_tc_model = stan_model("~/PhD/deconvolution/ARMET/inst/stan/ARMET_tc.stan")


	Sys.time() %>% print
browser()
	model  = switch(
		full_bayesian %>% `!` %>% sum(1),
		stanmodels$ARMET_tc,
		stanmodels$ARMET_tc_fix
	)

	fit =
		switch(
			full_bayes %>% `!` %>% sum(1),

			# HMC
			sampling(
				model, # ARMET_tc_model, #,
				chains=3, cores=3,
				iter=iterations, warmup=iterations-sampling_iterations,
				#control = list(stepsize = 0.05, adapt_delta = 0.99),
				#include = F, pars=c("prop_a", "prop_b", "prop_c", "prop_d", "prop_e"),
				#pars=c("prop_1", "prop_2", "prop_3", "prop_4", "exposure_rate", "lambda_log", "sigma_inv_log", "sigma_intercept_dec"),
				#,
				init = function () list(	lambda_log = lambda_log, sigma_inv_log = sigma_inv_log) # runif(G,  lambda_log - 1, lambda_log + 1)	)
				#save_warmup = FALSE,
				#pars = c("prop_1", "prop_2", "prop_3", "exposure_rate") #, "nb_sum") #,"mu_sum", "phi_sum"),
			) %>%
				{
					(.)  %>% rstan::summary() %$% summary %>% as_tibble(rownames="par") %>% arrange(Rhat %>% desc) %>% print
					(.)
				},

			vb_iterative(
				model, #ARMET_tc_model,
				output_samples=100,
				iter = 50000,
				tol_rel_obj=0.01,
				pars=c("prop_1", "prop_2", "prop_3","prop_4", "exposure_rate", "lambda_log", "sigma_inv_log", "sigma_intercept_dec"),
				#,
				init = function () list(	lambda_log = lambda_log) # runif(G,  lambda_log - 1, lambda_log + 1)	)

			)
		)

	Sys.time() %>% print

	# Produce results
	prop =
		fit %>%

		# If MCMC is used check divergencies as well
		ifelse_pipe(
			full_bayes,
			~ .x %>% parse_summary_check_divergence(),
			~ .x %>% parse_summary()
		) %>%

		# Parse
		separate(.variable, c(".variable", "level"), convert = T) %>%

		# Add tree information
		left_join(
			tree %>% data.tree::ToDataFrameTree("name", "C1", "C2", "C3", "C4") %>%
				as_tibble %>%
				select(-1) %>%
				rename(`Cell type category` = name) %>%
				gather(level, C, -`Cell type category`) %>%
				mutate(level = gsub("C", "", level)) %>%
				drop_na %>%
				mutate(C = C %>% as.integer, level = level %>% as.integer)
		) %>%

		# Add sample information
		left_join(
			df %>%
				filter(`query`) %>%
				distinct(Q, sample)
		)

	# Return
	list(

		# Matrix of proportions
		proportions =	prop,

		# Return the input itself
		input = input,

		# Return the fitted object
		fit = fit,

		# Return data source
		data_source = y_source,

		signatures = counts_baseline
	)

}



