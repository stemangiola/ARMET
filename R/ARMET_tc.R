

#' ARMET-tc main
#'
#' @description This function calls the stan model.
#'
#'
#' @importFrom tibble tibble
#'
#' @import dplyr
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
#' @param save_fit A boolean
#' @param seed An integer
#'
#' @return An ARMET object
#'
#' @export
#'
ARMET_tc = function(.data,
										.formula = ~ 1,
										.sample = NULL,
										.transcript = NULL,
										.abundance = NULL,
										reference = NULL,
										approximate_posterior = F,
										verbose =                           F,
										save_fit =                          F,
										seed =                              NULL,
										cores = 14,
										iterations = 250,
										sampling_iterations = 100,
										levels = 3,
										n_markers = my_n_markers ,
										do_regression = F) {

	# At the moment is not active
	full_bayesian = F

	input = c(as.list(environment()))

	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	col_names = get_sample_transcript_counts(.data, .sample, .transcript, .abundance)
	.sample = col_names$.sample
	.transcript = col_names$.transcript
	.abundance = col_names$.abundance


	# distinct_at is not released yet for dplyr, thus we have to use this trick
	df_for_edgeR <- .data %>%

		# Stop if any counts is NA
		error_if_counts_is_na(!!.abundance) %>%

		# Stop if there are duplicated transcripts
		error_if_duplicated_genes(!!.sample,!!.transcript,!!.abundance) %>%

		# Prepare the data frame
		select(!!.transcript,
					 !!.sample,
					 !!.abundance,
					 one_of(parse_formula(.formula))) %>%
		distinct() %>%

		# Check if data rectangular
		ifelse_pipe(
			(.) %>% check_if_data_rectangular(!!.sample,!!.transcript,!!.abundance, type = "soft") %>% `!` &
				TRUE, #!fill_missing_values,
			~ .x %>% eliminate_sparse_transcripts(!!.transcript)
		)

	# Create design matrix
	X =
		model.matrix(
			object = .formula,
			data = df_for_edgeR %>% select(!!.sample, one_of(parse_formula(.formula))) %>% distinct %>% arrange(!!.sample)
		)

	mix =
		.data %>%
		select(!!.sample,
					 !!.transcript,
					 !!.abundance,
					 one_of(parse_formula(.formula))) %>%
		distinct() %>%
		spread(!!.transcript, !!.abundance)

	shards = cores #* 2
	is_level_in = shards %>% `>` (0) %>% as.integer

	# Global properties - derived by previous analyses of the whole reference dataset
	sigma_intercept = 1.3420415
	sigma_slope = -0.3386389
	sigma_sigma = 1.1720851
	lambda_mu_mu = 5.612671
	lambda_sigma = 7.131593

	# Non centered
	lambda_mu_prior = c(6.2, 1)
	lambda_sigma_prior =  c(3.3 , 1)
	lambda_skew_prior =  c(-2.7, 1)
	sigma_intercept_prior = c(1.9 , 0.1)

	# Set up tree structure
	levels_in_the_tree = 1:4

	tree = 	data.tree::Clone(ARMET::tree) %>%	{
		# Filter selected levels
		data.tree::Prune(., function(x)
			x$level <= max(levels_in_the_tree) + 1)

		# Filter if not in referenc
		#data.tree::Prune(., function(x) ( x$name %in% (ARMET::ARMET_ref %>% distinct(`Cell type category`) %>% pull(1) %>% as.character) ))
		.
	}

	ct_in_nodes =
		tree %>%
		data.tree::ToDataFrameTree("name", "level", "C", "count", "isLeaf") %>%
		as_tibble %>%
		arrange(level, C) %>%
		filter(!isLeaf) %>%
		pull(count)

	# Get the number of leafs for every level
	ct_in_levels = foreach(l = levels_in_the_tree + 1, .combine = c) %do% {
		data.tree::Clone(tree) %>%
			ifelse_pipe((.) %>% data.tree::ToDataFrameTree("level") %>% pull(2) %>% max %>% `>` (l),
									~ {
										.x
										data.tree::Prune(.x, function(x)
											x$level <= l)
										.x
									})  %>%
			data.tree::Traverse(., filterFun = isLeaf) %>%
			length()
	}

	n_nodes = ct_in_nodes %>% length
	n_levels = ct_in_levels %>% length

	# Needed in the model
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

	Q = mix %>% nrow

	reference_filtered =
		ARMET::ARMET_ref %>%

		left_join(n_markers, by = c("ct1", "ct2")) %>%
		filter_reference(mix, n_markers) %>%
		select(-ct1,-ct2,-rank,-`n markers`) %>%
		distinct %>%

		# Select cell types in hierarchy
		inner_join(
			tree %>%
				data.tree::ToDataFrameTree("Cell type category", "C", "C1", "C2", "C3") %>%
				as_tibble %>%
				select(-1)

		)	%>%
		# Decrease the number of house keeping used
		anti_join({
			mdf = (.) %>%
				distinct(symbol, `house keeping`) %>%
				filter(`house keeping`)

			withr::with_seed(123, 	sample_frac(mdf, 0.5)) %>%
				distinct(symbol)
		})

	prop_posterior = get_null_prop_posterior(ct_in_nodes)

	######################################

	res1 = run_model(
		reference_filtered,
		mix,
		shards,
		1,
		full_bayesian,
		approximate_posterior,
		prop_posterior,
		iterations = iterations,
		sampling_iterations = sampling_iterations,
		X = X,
		do_regression = do_regression
	)

	df1 = res1[[1]]
	fit1 = res1[[2]]

	prop_posterior[[1]] = fit1 %>% draws_to_alphas("prop_1") %>% `[[` (1)

	draws_1 =
		fit1 %>%
		tidybayes::gather_draws(prop_1[Q, C1]) %>%
		ungroup() %>%
		select(-.variable) %>%
		mutate(.value_relative = .value)

	if (do_regression && length(parse_formula(.formula)) >0 ) {
		alpha_1 =
			fit1 %>%
			tidybayes::gather_draws(`alpha_[1]`[A, C], regex = T) %>%
			ungroup() %>%

			# rebuild the last component sum-to-zero
			rebuild_last_component_sum_to_zero %>%

			# Calculate relative 0 because of dirichlet relativity
			get_relative_zero %>%

			arrange(.chain, .iteration, .draw,     A) %>%

			nest(draws = -c(C, .variable, zero)) %>%

			left_join(
				tree %>%
					ToDataFrameTree("name", "level", sprintf("C%s", 1)) %>%
					filter(level == 1 + 1) %>%
					arrange(C1) %>%
					mutate(C = 1:n()) %>%
					select(name, C)  %>%
					rename(`Cell type category` = name)
			)
		alpha = alpha_1 %>% mutate(level = 1)
	}

	prop_1 =
		fit1 %>%
		tidybayes::gather_draws(`prop_[1]`[Q, C], regex = T) %>%
		drop_na  %>%
		ungroup() %>%

		# Add relative proportions
		mutate(.value_relative = .value) %>%

		# Add tree information
		left_join(
			tree %>% data.tree::ToDataFrameTree("name", "C1", "C2", "C3", "C4") %>%
				as_tibble %>%
				select(-1) %>%
				rename(`Cell type category` = name) %>%
				gather(level, C, -`Cell type category`) %>%
				mutate(level = gsub("C", "", level)) %>%
				filter(level == 1) %>%
				drop_na %>%
				mutate(C = C %>% as.integer, level = level %>% as.integer)
		) %>%

		# add sample annotation
		left_join(df1 %>% distinct(Q, sample), by = "Q")	%>%

		# If MCMC is used check divergencies as well
		ifelse_pipe(
			!approximate_posterior,
			~ .x %>% parse_summary_check_divergence(),
			~ .x %>% parse_summary() %>% rename(.value = mean)
		) %>%

		# Parse
		separate(.variable, c(".variable", "level"), convert = T) %>%



		# Add sample information
		left_join(df1 %>%
								filter(`query`) %>%
								distinct(Q, sample))

	prop = prop_1
	fit = list(fit1)
	df = list(df1)

	######################################

	if (levels > 1) {
		res2 = run_model(
			reference_filtered,
			mix,
			shards,
			2,
			full_bayesian,
			approximate_posterior,
			prop_posterior,
			draws_to_exposure(fit1)	,
			iterations = iterations,
			sampling_iterations = sampling_iterations,
			X = X,
			do_regression = do_regression
		)

		df2 = res2[[1]]
		fit2 = res2[[2]]

		prop_posterior[[2]] = fit2 %>% draws_to_alphas(sprintf("prop_%s", "a")) %>% `[[` (1)

		draws_a =
			fit2 %>%
			tidybayes::gather_draws(prop_a[Q, C]) %>%
			ungroup() %>%
			select(-.variable)

		draws_2 =
			draws_1 %>%
			left_join(
				fit2 %>%
					######## ALTERED WITH TREE
					tidybayes::gather_draws(`prop_[a]`[Q, C], regex = T) %>%
					###########################
				drop_na %>%
					ungroup() %>%
					left_join(
						######## ALTERED WITH TREE
						tibble(
							.variable = c("prop_a"),
							C1 = parents_lv2
						)
						##########################
					) %>%
						select(-.variable) %>%
							rename(.value2 = .value, C2 = C),
						by = c(".chain", ".iteration", ".draw", "Q",  "C1")
					) %>%
					group_by(.chain, .iteration, .draw, Q) %>%
					arrange(C1, C2) %>%
					mutate(
						C2 = tree$Get("C2") %>% na.omit,
						`Cell type category` = tree$Get("C2") %>% na.omit %>% names
					) %>%
					ungroup() %>%
					mutate(.value_relative = .value2) %>%
					mutate(.value2 = ifelse(.value2 %>% is.na, .value, .value * .value2))

				if (do_regression && length(parse_formula(.formula)) >0 ){
					alpha_2 =
						fit2 %>%
						tidybayes::gather_draws(`alpha_[2]`[A, C], regex = T) %>%
						ungroup() %>%

						# rebuild the last component sum-to-zero
						rebuild_last_component_sum_to_zero %>%

						# Calculate relative 0 because of dirichlet relativity
						get_relative_zero %>%

						arrange(.chain, .iteration, .draw,     A) %>%

						nest(draws = -c(C, .variable, zero)) %>%

						left_join(
							tree %>%
								ToDataFrameTree("name", "level", sprintf("C%s", 2)) %>%
								arrange(C2) %>%
								drop_na() %>%
								select(name, C2) %>%
								rename(`Cell type category` = name) %>%
								rename(C = C2)
						)
					alpha = alpha  %>% bind_rows(alpha_2 %>% mutate(level = 2))

				}

				prop_2 =
					draws_2 %>%
					select(.chain,
								 .iteration,
								 .draw,
								 Q,
								 C2 ,
								 `Cell type category`,
								 .value2,
								 .value_relative) %>%

					rename(C = C2, .value = .value2) %>%
					mutate(.variable = "prop_2") %>%
					mutate(level = 2) %>%

					# add sample annotation
					left_join(df2 %>% distinct(Q, sample), by = "Q")	%>%

					# If MCMC is used check divergencies as well
					ifelse_pipe(
						!approximate_posterior,
						~ .x %>% parse_summary_check_divergence(),
						~ .x %>% parse_summary() %>% rename(.value = mean)
					)

				# Eliminate
				prop = bind_rows(prop, prop_2)
				fit = fit %>% c(list(fit2))
				df = df %>% c(list(df2))
	}

	######################################

	if (levels > 2) {
		# browser()
		res3 = run_model(
			reference_filtered,
			mix,
			shards,
			3,
			full_bayesian,
			approximate_posterior,
			prop_posterior,
			draws_to_exposure(fit1),
			iterations = iterations,
			sampling_iterations = sampling_iterations	,
			X = X,
			do_regression = do_regression
		)

		df3 = res3[[1]]
		fit3 = res3[[2]]

		draws_3 =
			draws_2 %>%
			left_join(
				fit3 %>%
					######## ALTERED WITH TREE
					tidybayes::gather_draws(`prop_[b, c, d, e, f]`[Q, C], regex = T) %>%
					#########################
				drop_na %>%
					ungroup() %>%
					left_join(
						######## ALTERED WITH TREE
						tibble(
							.variable = c("prop_b", "prop_c", "prop_d", "prop_e", "prop_f"),
							C2 = parents_lv3
						)
						#########################
					) %>%
						select(-.variable) %>%
							rename(.value3 = .value, C3 = C),
						by = c(".chain", ".iteration", ".draw", "Q",  "C2")
					) %>%
					group_by(.chain, .iteration, .draw, Q) %>%
					arrange(C1, C2, C3) %>%
					mutate(
						C3 = tree$Get("C3") %>% na.omit,
						`Cell type category` = tree$Get("C3") %>% na.omit %>% names
					) %>%
					ungroup() %>%
					mutate(.value_relative = .value3) %>%
					mutate(.value3 = ifelse(.value3 %>% is.na, .value2, .value2 * .value3))

				if (do_regression && length(parse_formula(.formula)) >0 ){
					alpha_3 =
						fit3 %>%
						tidybayes::gather_draws(`alpha_[3]`[A, C], regex = T) %>%
						ungroup() %>%

						# rebuild the last component sum-to-zero
						rebuild_last_component_sum_to_zero %>%

						# Calculate relative 0 because of dirichlet relativity
						get_relative_zero %>%

						arrange(.chain, .iteration, .draw,     A) %>%

						nest(draws = -c(C, .variable, zero)) %>%

						left_join(
							tree %>%
								ToDataFrameTree("name", "level", sprintf("C%s", 3)) %>%
								arrange(C3) %>%
								drop_na() %>%
								select(name, C3) %>%
								rename(`Cell type category` = name) %>%
								rename(C = C3)
						)

					alpha = alpha  %>% bind_rows(alpha_3 %>% mutate(level = 3))

				}
				prop_3 =
					draws_3 %>%
					select(.chain,
								 .iteration,
								 .draw,
								 Q,
								 C3 ,
								 `Cell type category`,
								 .value3,
								 .value_relative) %>%
					rename(C = C3, .value = .value3) %>%
					mutate(.variable = "prop_3") %>%
					mutate(level = 3) %>%

					# add sample annotation
					left_join(df3 %>% distinct(Q, sample), by = "Q")	%>%

					# If MCMC is used check divergencies as well
					ifelse_pipe(
						!approximate_posterior,
						~ .x %>% parse_summary_check_divergence(),
						~ .x %>% parse_summary() %>% rename(.value = mean)
					) %>%

					left_join(df3 %>% distinct(Q, sample))

				prop = bind_rows(prop, prop_3)
				fit = fit %>% c(list(fit3))
				df = df %>% c(list(df3))

	}

	# Return
	list(
		# Matrix of proportions
		proportions =	prop %>%

			# Attach alpha if regression
			ifelse_pipe(
				do_regression && length(parse_formula(.formula)) >0 ,
				~ .x %>%
					nest(proportions = -c(`Cell type category`, C, level)) %>%
					left_join(
						alpha %>%	select(`Cell type category`, contains("alpha"), zero, level, draws),
						by = c("Cell type category", "level")
					)
			),

		# Return the input itself
		input = input,

		# Return the fitted object
		fit = fit,

		# # Return data source
		# data_source = y_source,
		signatures = df
	)

}


#' @export
run_model = function(reference_filtered,
										 mix,
										 shards,
										 lv,
										 full_bayesian,
										 approximate_posterior,
										 prop_posterior,
										 exposure_posterior = tibble(.mean = 0, .sd = 0)[0,],
										 iterations = 250,
										 sampling_iterations = 100,
										 X,
										 do_regression) {
	# Filter on level considered
	reference_filtered = reference_filtered %>% filter(level %in% lv)

	df = ref_mix_format(reference_filtered, mix)

	G = df %>% filter(!`query`) %>% distinct(G) %>% nrow()
	GM = df %>% filter(!`house keeping`) %>% distinct(symbol) %>% nrow()

	# For  reference MPI inference
	counts_baseline =
		df %>%

		# Eliminate the query part, not the house keeping of the query
		filter(!`query` | `house keeping`)  %>%

		format_for_MPI(shards)

	S = counts_baseline %>% distinct(sample) %>% nrow()
	N = counts_baseline %>% distinct(idx_MPI, `read count`, `read count MPI row`) %>%  count(idx_MPI) %>% summarise(max(n)) %>% pull(1)
	M = counts_baseline %>% distinct(start, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% max

	lambda_log = 	  counts_baseline %>% filter(!query) %>% distinct(G, lambda_log) %>% arrange(G) %>% pull(lambda_log)
	sigma_inv_log = counts_baseline %>% filter(!query) %>% distinct(G, sigma_inv_log) %>% arrange(G) %>% pull(sigma_inv_log)

	y_source =
		df %>%
		filter(`query` & !`house keeping`) %>%
		select(S, Q, `symbol`, `read count`, GM, sample) %>%
		left_join(
			df %>% filter(!query) %>% distinct(
				`symbol`,
				G,
				`Cell type category`,
				level,
				lambda_log,
				sigma_inv_log,
				GM,
				C
			)
		) %>%
		arrange(C, Q, symbol) %>%
		mutate(`Cell type category` = factor(`Cell type category`, unique(`Cell type category`)))

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
	S = counts_baseline_to_linear %>% distinct(S) %>% nrow

	# Counts idx for each level for each level
	counts_idx_lv_NA = counts_baseline_to_linear %>% filter(level %>% is.na) %>% pull(counts_idx)
	CL_NA = counts_idx_lv_NA %>% length

	# Level specific
	counts_idx_lv = counts_baseline_to_linear %>% filter(level == lv) %>% pull(counts_idx)
	CL_lv = counts_idx_lv %>% length

	# Deconvolution, get G only for markers of each level. Exclude house keeping
	G_lv_linear = counts_baseline %>% filter(level == lv) %>% select(G, GM, sprintf("C%s", lv)) %>% distinct() %>% arrange(GM,!!as.symbol(sprintf("C%s", lv))) %>% pull(G)
	G_lv = G_lv_linear %>% length

	# Observed mix counts
	y_linear_lv = y_source %>% filter(level == lv) %>% distinct(GM, Q, S, `read count`) %>% arrange(GM, Q) %>% pull(`read count`)

	# Observed mix samples indexes
	y_linear_S_lv = y_source %>% filter(level == lv) %>% distinct(GM, Q, S, `read count`) %>% arrange(GM, Q) %>% pull(S)

	# Lengths indexes
	Y_lv = y_linear_lv %>% length


	MPI_data = get_MPI_df(counts_baseline_to_linear,
												y_source,
												counts_baseline,
												shards,
												lv)

	# Dirichlet regression
	A = X %>% ncol



	model  = switch(full_bayesian %>% `!` %>% sum(1),
									stanmodels$ARMET_tc,
									stanmodels$ARMET_tc_fix)

	additional_par_to_save  = switch(full_bayesian %>% `!` %>% sum(1),
																	 c("lambda_log", "sigma_inv_log"),
																	 c())

	# library(rstan)
	# fileConn<-file("~/.R/Makevars")
	# writeLines(c( "CXX14FLAGS += -O2","CXX14FLAGS += -DSTAN_THREADS", "CXX14FLAGS += -pthread"), fileConn)
	# close(fileConn)
	# ARMET_tc_model = rstan::stan_model("~/PhD/deconvolution/ARMET/inst/stan/ARMET_tc_fix.stan", auto_write = F)

	exposure_rate_init = switch(
		(lv > 1) %>% `!` %>% as.numeric %>% sum(1),
		exposure_posterior %>% pull(1),
		runif(S, -0.5, 0.5)
	) %>% as.array

	exposure_rate_multiplier = sd(exposure_rate_init) %>% ifelse_pipe((.) %>% is.na, ~ 0.1)

	init_list = list(lambda_log = lambda_log,
									 sigma_inv_log = sigma_inv_log) %>%
		ifelse_pipe(!full_bayesian,
								~ .x %>% c(list(exposure_rate = exposure_rate_init)))


	Sys.setenv("STAN_NUM_THREADS" = shards)

	#if(lv == 3) browser()

	list(df,
			 switch(
			 	approximate_posterior %>% sum(1),

			 	# HMC
			 	sampling(
			 		model,
			 		#ARMET_tc_model, #,
			 		chains = 3,
			 		cores = 3,
			 		iter = iterations,
			 		warmup = iterations - sampling_iterations,
			 		data = MPI_data %>% c(prop_posterior),
			 		# pars=
			 		# 	c("prop_1", "prop_2", "prop_3", sprintf("prop_%s", letters[1:9])) %>%
			 		# 	c("alpha_1", sprintf("alpha_%s", letters[1:9])) %>%
			 		# 	c("exposure_rate") %>%
			 		# 	c("lambda_UFO") %>%
			 		# 	c("prop_UFO") %>%
			 		# 	c(additional_par_to_save),
			 		init = function ()
			 			init_list,
			 		save_warmup = FALSE
			 	) %>%
			 		{
			 			(.)  %>% rstan::summary() %$% summary %>% as_tibble(rownames = "par") %>% arrange(Rhat %>% desc) %>% print
			 			(.)
			 		},

			 	vb_iterative(
			 		model,
			 		#ARMET_tc_model,
			 		output_samples = 100,
			 		iter = 50000,
			 		tol_rel_obj = 0.01,
			 		data = MPI_data,
			 		pars = c(
			 			"prop_1",
			 			"prop_2",
			 			"prop_3",
			 			"prop_4",
			 			"exposure_rate",
			 			"lambda_log",
			 			"sigma_inv_log",
			 			"sigma_intercept_dec"
			 		),
			 		#,
			 		init = function ()
			 			list(lambda_log = lambda_log, sigma_inv_log = sigma_inv_log) # runif(G,  lambda_log - 1, lambda_log + 1)	)

			 	)
			 ))

}

#' @export
test_differential_composition = function(.data, credible_interval = 0.95) {

	.data$proportions %>%
		mutate(regression = map2(draws, zero,
			~ .x %>%
				group_by(A) %>%
				tidybayes::median_qi(.width = credible_interval) %>%
				ungroup()  %>%
				pivot_wider(
					names_from = A,
					values_from = c(.value, .lower, .upper),
					names_prefix = "alpha"
				) %>%
				mutate(
					.value_alpha2 = .value_alpha2 - .y,
					.lower_alpha2 = .lower_alpha2 - .y,
					.upper_alpha2 = .upper_alpha2 - .y
				) %>%
				mutate(significant = (.lower_alpha2 * .upper_alpha2) > 0)
		)) %>%
		unnest(cols = c(regression))

}

#' Create polar plot of results
#' @rdname ARMET_plotTree
#'
#' Prints a report of the hipothesis testing
#'
#' @import ggplot2
#' @import tibble
#' @import dplyr
#' @import data.tree
#' @importFrom ape as.phylo
#'
#'
#' @param ARMET-tc object
#'
#' @return a ggplot
#'
#' @export
plot_polar = function(	.data,
												size_geom_text = 3.5,
												my_breaks=c(0, 0.01, 0.03, 0.05, 0.1, 0.3, 0.5, 0.7, 1),
												prop_filter = 0.005,
												barwidth = 0.5,
												barheight = 2,
												legend_justification = 0.67,
												fill_direction = 1){


	# Create annotation
	internal_branch_length = 40
	external_branch_length = 10


	dd =

		# Integrate data
		ARMET::tree %>%
		data.tree::ToDataFrameTree("name", "isLeaf", "level", "leafCount", pruneFun = function(x)	x$level <= 4) %>%
		as_tibble() %>%
		rename(`Cell type category` = name) %>%
		mutate(level = level -1) %>%
		left_join(.data) %>%

		# process
		mutate(Estimate = ifelse(significant, .value_alpha2, NA)) %>%
		mutate(branch_length = ifelse(isLeaf, 0.1, 2)) %>%

		# Correct branch length
		mutate(branch_length = ifelse(!isLeaf, internal_branch_length,	external_branch_length) ) %>%
		mutate(branch_length =  ifelse(	isLeaf, branch_length + ((max(level) - level) * internal_branch_length),	branch_length	))



	xx = dd %>%
		mutate(leafCount_norm = leafCount/max(leafCount)) %>%
		group_by(level) %>%
		arrange(desc(leafCount>1)) %>%
		mutate(leafCount_norm_cum = cumsum(leafCount_norm)) %>%
		ungroup() %>%
		arrange(level) %>%
		mutate(length_error_bar = leafCount_norm - 0.005) %>%
		mutate(x = leafCount_norm_cum - 0.5 * leafCount_norm) %>%
		mutate(angle = x * -360) %>%
		mutate(angle = ifelse(angle < - 240 | angle > -120, angle, angle - 180)) %>%
		mutate(angle = ifelse(angle > - 360, angle, angle + 360)) %>%
		#mutate(median_proportion = ifelse(level<max(level), median_proportion + (0.008*level), median_proportion  )) %>%
		filter(level > 0) %>%

		# Calculate median proportions
		filter(.value_alpha1 %>% is.na %>% `!`) %>%
		mutate(summarised_proportion =
					 	proportions %>% map(~ .x %>% summarise(median_proportion = .value %>% median))
		) %>%
		unnest(summarised_proportion)
		#mutate(median_proportion = median_proportion / max(median_proportion))

	my_rescale = function(x) { as.character( format(round( x * xx %>% pull(median_proportion) %>% max, 2), nsmall = 2)) }
	cusotm_root_trans = function() scales::trans_new("cusotm_root",function(x) x^(1/4), function(x) x^(4))
	DTC_scale =  xx %>% pull(Estimate) %>% abs() %>% { switch( (!all(is.na(.))) + 1, 1, .) } %>% max(na.rm = T)
	y_text_out = 5.0625
	y_text_in = 1.4641


	xx %>%
		{
			# Case if none is significant
			switch(
				(!length(na.omit(xx$Estimate))>0) + 1,

				# If there are significant
				ggplot(data=(.), aes(x = x,fill = Estimate,size = 1/sqrt(level))),

				# If there are not
				ggplot(data=(.),aes(x = x,fill = "Non significant",size = 1/sqrt(level)))
			)
		}	+
		annotate("rect", xmin=0, xmax=1, ymin=1, ymax= 2.0736, fill="grey95") +
		annotate("rect", xmin=0, xmax=1, ymin=2.0736, ymax=3.8416, fill="grey90") +
		annotate("rect", xmin=0, xmax=1, ymin=3.8416, ymax=6.5536, fill="grey87") +
		geom_bar(
			data = xx %>% filter(level == 1),
			aes(width = leafCount_norm, y = `median_proportion`), # / leafCount_norm),
			color = "grey20",  stat = "identity"
		) +
		geom_errorbar(
			data = xx %>% filter(level == 1) %>% mutate(x = ifelse(`Cell type category`=="immune_cell", 0.5, x)),
			aes(width = 0 , ymin=`median_proportion`, ymax=y_text_out-0.7), # / leafCount_norm),
			color = "grey20",  stat = "identity"
		) +
		geom_text(
			data =
				xx %>%
				filter(level == 1) %>%
				mutate(x = ifelse(`Cell type category`=="immune_cell", 0.5, x)) %>%
				mutate(angle = ifelse(`Cell type category`=="immune_cell", -0, angle)) %>%
				mutate(`Cell type category` = sub("\\s+$", "", `Cell type category`)),
			aes(label=`Cell type category`, y = y_text_out, angle= angle  ) ,size =size_geom_text ) +
		#scale_x_continuous(labels = xx %>% filter(level == 1) %>% pull(`Cell type category`), breaks = xx %>% filter(level == 1) %>% pull(leafCount_norm_cum) - 0.5 * xx %>% filter(level == 1) %>% pull(leafCount_norm)) +
		geom_bar(
			data = xx %>% filter(level == 2),
			aes(width = leafCount_norm, y = `median_proportion` ), # / leafCount_norm),
			color = "grey20",  stat = "identity"
		) +
		geom_errorbar(
			data = xx %>% filter(level == 2),
			aes(width = 0 , ymin=`median_proportion`, ymax=2.2736), # / leafCount_norm),
			color = "grey20",  stat = "identity",linetype="dotted"
		) +
		geom_text(
			data =
				xx %>%
				filter(level == 2) %>%
				mutate(`Cell type category` = sub("\\s+$", "", `Cell type category`)),
			aes(label=`Cell type category`, y = 2.8561, angle= angle) ,size =size_geom_text) +
		#scale_x_continuous(labels = xx %>% filter(level == 2) %>% pull(`Cell type category`), breaks = xx %>% filter(level == 2) %>% pull(leafCount_norm_cum) - 0.5 * xx %>% filter(level == 2) %>% pull(leafCount_norm)) +
		geom_bar(
			data = xx %>% filter(level == 3),
			aes(width = leafCount_norm, y = `median_proportion` ), # / leafCount_norm),
			color = "grey20", stat = "identity"
		)  +
		{
			# Make plotting robust if no level 3 cell types were detected
			switch(
				(!  xx %>% filter(level == 3) %>% filter(`median_proportion`>prop_filter | !is.na(Estimate)) %>% nrow() > 0) + 1,

				# If there are cell types
				geom_errorbar(
					data = xx %>% filter(level == 3) %>% filter(`median_proportion`>prop_filter | !is.na(Estimate)),
					aes(width = 0 , ymin=`median_proportion`, ymax=y_text_in-0.2), # / leafCount_norm),
					color = "grey20",  stat = "identity",linetype="dotted"
				) ,

				# If there are NOT cell types
				geom_errorbar(ymin=0, ymax=0)
			)
		} +
		{
			# Make plotting robust if no level 3 cell types were detected
			switch(
				(!  xx %>% filter(level == 3) %>% filter(`median_proportion`>prop_filter | !is.na(Estimate)) %>% nrow() > 0) + 1,

				# If there are cell types

				geom_text(
					data =
						xx %>%
						filter(level == 3)  %>%
						filter(
							`median_proportion`>prop_filter |
								!is.na(Estimate)
						) %>%
						mutate(`Cell type category` = sub("\\s+$", "", `Cell type category`)),
					aes(label=`Cell type category`, y = y_text_in , angle= angle) ,size =size_geom_text
				),

				# If there are NOT cell types
				geom_errorbar(ymin=0, ymax=0)
			)
		} +
		{
			# Case if none is significant
			switch(
				(!length(na.omit(xx$Estimate))>0) + 1,

				# If there are significant
				scale_fill_distiller(
					palette = "Spectral",
					na.value = 'white',
					direction = fill_direction, name = "Trend",
					limits=c(-DTC_scale, DTC_scale)
				) ,

				# If there are not
				scale_fill_manual(values= c("Non significant" = "white" ), guide=FALSE)
			)

		} +
		{
			# Case if none is significant
			switch(
				(!length(na.omit(xx$Estimate))>0) + 1,

				# If there are significant
				guides(
					fill = guide_colorbar(
						label.position = "left",
						title.position = "left",
						title.hjust = 0.5,
						override.aes=list(fill=NA),
						ticks.colour = "black",
						barwidth = barwidth,
						barheight = barheight
					)
				) ,

				# If there are not
				scale_fill_manual(values= c("Non significant" = "white" ), guide=FALSE)
			)

		} +
		scale_y_continuous( breaks=my_breaks, trans = cusotm_root_trans(), labels = my_rescale ) +
		scale_size(range = c(0.3, 0.8), guide=FALSE) +
		theme_bw() +
		theme(
			axis.text.x = element_blank(),
			axis.ticks.x.top = element_blank(),
			axis.title.x=element_blank(),
			axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0), hjust = 0.65, vjust = 1),
			panel.border = element_blank(),
			axis.line.y = element_line(),
			panel.grid  = element_blank(),
			legend.position=c(0,0.05),
			legend.justification=c(legend_justification, 0),
			legend.title=element_text(angle = 90),
			legend.background = element_rect(colour = "transparent", fill = alpha("red", 0))
		)	+ ylab("Cell type proportion") +
		coord_polar(theta = "x")
}


