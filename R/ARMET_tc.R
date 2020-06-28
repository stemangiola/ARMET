#' ARMET_tc_continue
#' 
#' @description This function
#' 
#' @export
ARMET_tc_continue = function(armet_obj, level, model = stanmodels$ARMET_tc_fix_hierarchical){
	
	internals = armet_obj$internals
	input = armet_obj$input
	
	
	res = run_model(
		reference_filtered = internals$reference_filtered,
		mix = internals$mix,
		shards = input$cores,
		lv = level,
		full_bayesian = input$full_bayesian,
		approximate_posterior = input$approximate_posterior,
		internals$prop_posterior,
		exposure_posterior = when(
			level,
			(.)==1 ~ tibble(.mean = 0, .sd = 0)[0,], 
			~ draws_to_exposure(internals$fit[[1]])
		),
		iterations = input$iterations,
		sampling_iterations = input$sampling_iterations	,
		X = internals$X,
		do_regression = input$do_regression,
		family = input$family,
		cens = internals$cens,
		tree_properties = internals$tree_properties,
		Q = internals$Q,
		model = model,
		prior_survival_time = internals$prior_survival_time
	)
	
	df = res[[1]]
	fit = res[[2]]
	
	fit_prop_parsed = 
		fit %>%
		draws_to_tibble("prop_", "Q", "C") %>%
		filter(!grepl("_UFO|_rng", .variable))  %>%
		mutate(Q = Q %>% as.integer)
	
	draws = get_draws(fit_prop_parsed, level, internals)	
	
	prop = get_props(draws, level, df, input$approximate_posterior)	
	
	internals$prop = bind_rows(internals$prop , prop) 
	internals$fit = internals$fit %>% c(list(fit))
	internals$df = internals$df %>% c(list(df))
	internals$prop_posterior[sprintf("%s_prior", fit_prop_parsed %>% distinct(.variable) %>% pull())] = fit_prop_parsed %>% group_by(Q, C, .variable) %>% prop_to_list
	internals$draws = internals$draws %>% c(list(draws))
	
	if (input$do_regression && length(parse_formula(input$.formula)$covariates) >0 )
		internals$alpha = internals$alpha  %>% bind_rows( get_alpha(fit, level, input$family) )
	 
	# Return
	list(
		# Matrix of proportions
		proportions =	
			input$.data %>%
			select(c(sample, (.) %>% get_specific_annotation_columns(sample))) %>%
			distinct() %>%
			left_join(internals$prop) %>%
			
			# Attach alpha if regression
			ifelse_pipe(
				input$do_regression && length(parse_formula(input$.formula)$covariates) >0 ,
				~ .x %>%
					nest(proportions = -c(`Cell type category`, C, level)) %>%
					left_join(
						internals$alpha %>%	select(`Cell type category`, contains("alpha"), level, draws, rng, .variable, Rhat),
						by = c("Cell type category", "level")
					)
			),
		
		# # Return the input itself
		input = input,
		
		# Return the fitted object
		internals = internals
	)
	
}

#' ARMET-tc main
#'
#' @description This function calls the stan model.
#'
#'
#'
#' @import dplyr
#' @import purrr
#'
#' @importFrom tidyr spread
#' @importFrom tidyr gather
#' @importFrom tidyr drop_na
#' @importFrom rlang enquo
#' @importFrom tibble tibble
#' @importFrom stringr str_split
#' @importFrom tidybayes gather_draws
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
										family = "dirichlet",
										reference = NULL,
										approximate_posterior = F,
										verbose =                           F,
										save_fit =                          F,
										seed =                              NULL,
										cores = 14,
										iterations = 250,
										sampling_iterations = 100,
										levels = 1,
										.n_markers = n_markers ,
										do_regression = T, 
										prior_survival_time = c(),
										model = stanmodels$ARMET_tc_fix_hierarchical) {

	# At the moment is not active
	full_bayesian = F

	input = c(as.list(environment()))
	input$.formula = .formula
	
	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	col_names = get_sample_transcript_counts(.data, .sample, .transcript, .abundance)
	.sample = col_names$.sample
	.transcript = col_names$.transcript
	.abundance = col_names$.abundance

	# Rename columns mix
	.data = .data %>% rename( sample = !!.sample, symbol = !!.transcript ,  count = !!.abundance)
	input$.data = .data
	
	
	# Warning is sensitive names in columns
	names_taken = c("level") 
	if(.data %>% colnames %in% names_taken %>% any) stop(sprintf("ARMET says: your input data frame includes reserved column names: %s", names_taken))
	
	# Checkif count is integer
	if(.data %>% select(count) %>% lapply(class) %>% unlist() %>% equals("integer") %>% `!`)
		stop(sprintf("ARMET says: the %s column must be integer as the deconvolution model is Negative Binomial", quo_name(.abundance)))

	# Check family
	#if(family %in% c("dirichlet", "beta") %>% any %>% `!`) stop("ARMET says: Please choose between dirichlet or beta families")

	# Covariate column
	cov_columns =
		parse_formula(.formula)$covariates %>%
		map_chr(~ .x %>% gsub("censored\\(|\\)| ", "", .) %>% str_split("\\,") %>% `[[` (1) %>% `[` (1)) %>%
		ifelse_pipe((.) %>% is.null, ~ c())

	# Do regresson
	if(length(cov_columns) > 0 & (cov_columns %>% is.na %>% `!`)) do_regression = T

	# distinct_at is not released yet for dplyr, thus we have to use this trick
	df_for_edgeR <- .data %>%

		# Stop if any counts is NA
		error_if_counts_is_na(count) %>%

		# Stop if there are duplicated transcripts
		error_if_duplicated_genes(sample,symbol,count) %>%

		# Prepare the data frame
		select(symbol,
					 sample,
					 count,
					 one_of(cov_columns)) %>%
		distinct() %>%

		# Check if data rectangular
		ifelse_pipe(
			(.) %>% check_if_data_rectangular(sample,symbol,count, type = "soft") %>% `!` &
				TRUE, #!fill_missing_values,
			~ .x %>% eliminate_sparse_transcripts(symbol)
		)

	# Censoring column
	.cens_column = parse_formula(.formula)$covariates %>% grep("censored(", ., fixed = T, value = T)  %>% gsub("censored\\(|\\)| ", "", .) %>% str_split("\\,") %>% ifelse_pipe(length(.)>0, ~.x %>% `[[` (1) %>% `[` (-1), ~NULL)
	.cens_value_column = parse_formula(.formula)$covariates %>% grep("censored(", ., fixed = T, value = T)  %>% gsub("censored\\(|\\)| ", "", .) %>% str_split("\\,") %>% ifelse_pipe(length(.)>0, ~.x %>% `[[` (1) %>% `[` (1), ~NULL)
	
	if(length(.cens_column) == 1) {
		cens = .data %>% select(sample, .cens_column) %>% distinct %>% arrange(sample) %>% pull(2)
		
		# Check cens right type
		if(typeof(cens) %in% c("integer", "logical") %>% any %>% `!`) stop("ARMET says: censoring variable should be logical of integer (0,1)")
		if(length(prior_survival_time) == 0) stop("AMET says: you really need to provide third party survival time for your condition/disease")
		
		sd_survival_months = .data %>%  select(sample, .cens_value_column) %>% distinct %>% pull(.cens_value_column) %>% sd
		df_for_edgeR = df_for_edgeR %>% mutate(!!.cens_value_column := !!as.symbol(.cens_value_column) / sd_survival_months)
		prior_survival_time = prior_survival_time / sd_survival_months

	}
	else{
		cens = NULL
	} 
	
	# Create design matrix
	if(length(cov_columns) > 0 & (cov_columns %>% is.na %>% `!`)) my_formula = as.formula( paste("~",paste(cov_columns, collapse = "+")))
	else my_formula = .formula

	X =
		model.matrix(
			object = 	my_formula,
			data = df_for_edgeR %>% select(sample, one_of(cov_columns)) %>% distinct %>% arrange(sample)
		)

	mix =
		.data %>%
		select(sample, symbol, count, one_of(parse_formula(.formula)$covariates)) %>%
		distinct() %>%
		spread(symbol, count)

	shards = cores #* 2
	is_level_in = shards %>% `>` (0) %>% as.integer

	tree = 	data.tree::Clone(ARMET::tree) 

	
	# Print overlap descriptive stats
	#get_overlap_descriptive_stats(mix %>% slice(1) %>% gather(symbol, count, -sample), reference)

	# Prepare data frames -
	# For Q query first
	# For G house keeing first
	# For GM level 1 first

	Q = mix %>% nrow

	reference_filtered =
		ARMET::ARMET_ref %>%

		left_join(.n_markers, by = c("ct1", "ct2")) %>%
		filter_reference(mix, .n_markers) %>%
		select(-ct1,-ct2,-rank,-`n markers`) %>%
		distinct %>%

		# Select cell types in hierarchy
		inner_join(
			tree %>%
				data.tree::ToDataFrameTree("Cell type category", "C", "C1", "C2", "C3", "C4") %>%
				as_tibble %>%
				select(-1),
			by = "Cell type category"
		)	%>%
		
		# Decrease the number of house keeping used
		anti_join({
			mdf = (.) %>%
				distinct(symbol, `house keeping`) %>%
				filter(`house keeping`)

			withr::with_seed(123, 	sample_frac(mdf, 0.5)) %>%
				distinct(symbol)
		},
		by = "symbol"
	)

	tree_propeties = get_tree_properties(tree)
	
	# Default internals
	internals = 
		list(
			prop = NULL,
			fit = NULL,
			df = NULL,
			prop_posterior = get_null_prop_posterior(tree_propeties$ct_in_nodes),
			alpha = NULL,
			Q = Q,
			reference_filtered = reference_filtered,
			mix = mix,
			X = X,
			cens = cens,
			tree_properties = tree_propeties,
			prior_survival_time = prior_survival_time
		) 
	
	internals = 
		run_lv_1(
			internals,
			shards,
			levels,
			full_bayesian,
			approximate_posterior,
			iterations = iterations,
			sampling_iterations = sampling_iterations	,
			do_regression = do_regression,
			family = family,
			.formula = .formula,
			model = model
		)
	 
	# Return
	list(
		# Matrix of proportions
		proportions =	
			.data %>%
			select(c(sample, (.) %>% get_specific_annotation_columns(sample))) %>%
			distinct() %>%
			left_join(internals$prop) %>%

			# Attach alpha if regression
			ifelse_pipe(
				do_regression && length(parse_formula(.formula)$covariates) >0 ,
				~ .x %>%
					nest(proportions = -c(`Cell type category`, C, level)) %>%
					left_join(
						internals$alpha %>%	select(`Cell type category`, contains("alpha"), level, draws, rng, .variable, Rhat),
						by = c("Cell type category", "level")
					)
			),

		# # Return the input itself
		input = input,

		# Return the fitted object
		internals = internals
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
										 do_regression,
										 family = "dirichlet",
										 cens,
										 tree_properties,
										 Q,
										 model = stanmodels$ARMET_tc_fix_hierarchical,
										 prior_survival_time = c()) {
	
	Q = Q

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
	N = counts_baseline %>% distinct(idx_MPI, count, `read count MPI row`) %>%  count(idx_MPI) %>% summarise(max(n)) %>% pull(1)
	M = counts_baseline %>% distinct(start, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% max

	lambda_log = 	  counts_baseline %>% filter(!query) %>% distinct(G, lambda_log) %>% arrange(G) %>% pull(lambda_log)
	sigma_inv_log = counts_baseline %>% filter(!query) %>% distinct(G, sigma_inv_log) %>% arrange(G) %>% pull(sigma_inv_log)

	y_source =
		df %>%
		filter(`query` & !`house keeping`) %>%
		select(S, Q, symbol, count, GM, sample) %>%
		left_join(
			df %>% filter(!query) %>% distinct(
				symbol,
				G,
				`Cell type category`,
				level,
				lambda_log,
				sigma_inv_log,
				GM,
				C
			),
			by = c("symbol", "GM")
		) %>%
		arrange(C, Q, symbol) %>%
		mutate(`Cell type category` = factor(`Cell type category`, unique(`Cell type category`)))

	counts_baseline_to_linear =
		counts_baseline %>%
		filter_house_keeping_query_if_fixed(full_bayesian) %>%
		arrange(G, S) %>%
		mutate(counts_idx = 1:n()) %>%
		mutate(S = S %>% as.factor %>% as.integer)

	counts_linear = counts_baseline_to_linear %>%  pull(count)
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
	y_linear_lv = y_source %>% filter(level == lv) %>% distinct(GM, Q, S, count) %>% arrange(GM, Q) %>% pull(count)

	# Observed mix samples indexes
	y_linear_S_lv = y_source %>% filter(level == lv) %>% distinct(GM, Q, S, count) %>% arrange(GM, Q) %>% pull(S)

	# Lengths indexes
	Y_lv = y_linear_lv %>% length


	MPI_data = get_MPI_df(counts_baseline_to_linear,
												y_source,
												counts_baseline,
												shards,
												lv)

	# Dirichlet regression
	A = X %>% ncol



	
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

	fam_dirichlet = family == "dirichlet"

	if(cens %>% is.null) cens =  rep(0, Q)
	which_cens = which(cens == 1)  %>% as.array()
	which_not_cens = which(cens == 0) %>% as.array()
	how_many_cens = length(which_cens)

	max_unseen = ifelse(how_many_cens>0, max(X[,2]), 0 )
	if(is.null(prior_survival_time)) prior_survival_time = array(1)[0]
	spt = length(prior_survival_time)

	
	# model  = stanmodels$ARMET_tc_fix_hierarchical
	# switch(fam_dirichlet %>% `!` %>% sum(1),
	# 								stanmodels$ARMET_tc_fix_hierarchical,
	# 								stanmodels$ARMET_tc_fix)

	fit = 
		sampling(
			model,
			#rstan::stan_model("~/PhD/deconvolution/ARMET/inst/stan/ARMET_tc_fix_hierarchical.stan", auto_write = F),
			chains = 3,
			cores = 3,
			iter = iterations,
			warmup = iterations - sampling_iterations,
			data = MPI_data %>% c(prop_posterior) %>% c(tree_properties),
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
			(.)  %>% rstan::summary() %$% summary %>% as_tibble(rownames = "par") %>% arrange(Rhat %>% desc) %>% filter(Rhat > 1.5) %>% ifelse_pipe(nrow(.) > 0, ~ .x %>% print)
			(.)
		}

	list(df,
			 switch(
			 	approximate_posterior %>% sum(1),

			 	fit,

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
get_signatures = function(.data){
	.data$proportions %>%
		filter(.variable %>% is.na %>% `!`) %>%
		select(-proportions, -rng) %>%
		unnest(draws) %>%
		
		# Group
		nest(node = -c(level, .variable, A)) %>%
		mutate(node = map(
			node,
			~ .x %>%

				# Build combination of cell types
				nanny::permute_nest(
					.names_from = `Cell type category`,
					.values_from = c(.value)
				) %>% 
				
				# Perform calculation
				mutate(prob_df = map(
					data, 
					~.x %>%
						nest(data = -c(`Cell type category`)) %>% 
						mutate(med = map_dbl(data, ~.x$.value %>% median )) %>% 
						mutate(med = rev(med)) %>% 
						mutate(frac_up = map2_dbl(data, med, ~ (.x$.value  > .y) %>% sum %>% magrittr::divide_by(nrow(.x)))) %>% 
						mutate(frac_down = map2_dbl(data, med, ~ (.x$.value  < .y) %>% sum %>% magrittr::divide_by(nrow(.x)))) %>%
						mutate(prob = ifelse(frac_up > frac_down, frac_up, frac_down)) %>%
						
						# stretch to o 1 interval
						mutate(prob = (prob-0.5)/0.5) %>%
						
						# Insert sign
						mutate(prob = ifelse(frac_up < frac_down, -prob, prob)) %>%
						select(`Cell type category`, prob)
				)) %>%
				select(-data) %>%
				unnest(prob_df) %>% 
				filter(`Cell type category_1` == `Cell type category`) %>%
				select(-`Cell type category`)
				
		)) %>%
		unnest( node)
}

#' @export
test_differential_composition = function(.data, credible_interval = 0.90, cluster_CI = 0.55) {

 
	.d = 
		.data$proportions %>%
		filter(.variable %>% is.na %>% `!`) %>%
		cluster_posterior_slopes(credible_interval = cluster_CI) %>%
		extract_CI(credible_interval)
	
	dx = list()	
	
	# Level 1
	if(.d %>% filter(level ==1) %>% nrow %>% `>` (0))
	dx = 
		.d %>%
		
		filter(level ==1) %>%
		mutate(fold_change_ancestor = 0) %>%
		identify_baseline_by_clustering( ) %>%
		
		mutate(significant = ((.lower_alpha2 - zero) * (.upper_alpha2 - zero)) > 0) %>%
		mutate(fold_change  = ifelse(significant, .value_alpha2, 0))
		
	# Level 2
	if(.d %>% filter(level ==2) %>% nrow %>% `>` (0))
		dx =	dx %>% bind_rows(
				.d %>%
		
		filter(level ==2) %>%
		left_join(ancestor_child, by = "Cell type category") %>%
		left_join( dx %>%	select(ancestor = `Cell type category`, fold_change_ancestor = fold_change),  by = "ancestor" )  %>%
		identify_baseline_by_clustering( ) %>%
		
		mutate(significant = ((.lower_alpha2 - zero) * (.upper_alpha2 - zero)) > 0) %>%
		mutate(fold_change  = ifelse(significant, .value_alpha2, 0))
		) 
	
	
	# Level 3
	if(.d %>% filter(level ==3) %>% nrow %>% `>` (0))
		dx =	dx %>% bind_rows(
		.d %>%
		
		filter(level ==3) %>%
		left_join(ancestor_child, by = "Cell type category") %>%
		left_join( dx %>%	select(ancestor = `Cell type category`, fold_change_ancestor = fold_change) ,  by = "ancestor")  %>%
		identify_baseline_by_clustering( ) %>%
		
		mutate(significant = ((.lower_alpha2 - zero) * (.upper_alpha2 - zero)) > 0) %>%
		mutate(fold_change  = ifelse(significant, .value_alpha2, 0))
		)
	
	# Level 4
	if(.d %>% filter(level ==4) %>% nrow %>% `>` (0))
		dx =	dx %>% bind_rows(
			.d %>%
				
				filter(level ==4) %>%
				left_join(ancestor_child, by = "Cell type category") %>%
				left_join( dx %>%	select(ancestor = `Cell type category`, fold_change_ancestor = fold_change) ,  by = "ancestor")  %>%
				identify_baseline_by_clustering( ) %>%
				
				mutate(significant = ((.lower_alpha2 - zero) * (.upper_alpha2 - zero)) > 0) %>%
				mutate(fold_change  = ifelse(significant, .value_alpha2, 0))
		)
	
	dx
	
	
	# # Get zero
	# .data$proportions %>%
	# 	filter(.variable %>% is.na %>% `!`) %>%
	# 	extract_CI %>%
	# 	
	# 	# Link zero to closest posterior
	# 	nest(data = -.variable) %>%
	# 	mutate(
	# 		data = map(
	# 			data, 
	# 			~ .x %>%
	# 				mutate(
	# 					zero = 
	# 						.x %>%
	# 						mutate(diff = abs(zero - .value_alpha2)) %>%
	# 						arrange(diff) %>% 
	# 						slice(1) %>%
	# 						pull(.value_alpha2)
	# 				)
	# 			)
	# 	) %>%
	# 	unnest(data) %>%
	# 	
	# 	# Label signifcant
	# 	# mutate(
	# 	# 	.value_alpha2 = .value_alpha2 - zero,
	# 	# 	.lower_alpha2 = .lower_alpha2 - zero,
	# 	# 	.upper_alpha2 = .upper_alpha2 - zero
	# 	# ) %>%
	# 	mutate(significant = ((.lower_alpha2 - zero) * (.upper_alpha2 - zero)) > 0) %>%
	# 	
	# 	# Adjust ifancestor is significant
	# 	# equential important becaue we hav toudatecildren on level at the time
	# 	left_join(ancestor_child)	%>%
	# 	left_join( (.) %>% distinct(ancestor = `Cell type category`, fold_change_ancestor = ifelse(significant, .value_alpha2, 0)) )	%>%
	# 	group_by(.variable) %>%
	# 	mutate(zero = ifelse(level == 2 & fold_change_ancestor > 0, min(.value_alpha2), zero)) %>%
	# 	mutate(zero = ifelse(level == 2 & fold_change_ancestor < 0, max(.value_alpha2), zero)) %>%
	# 	mutate(zero = ifelse(level == 3 & fold_change_ancestor > 0, min(.value_alpha2), zero)) %>%
	# 	mutate(zero = ifelse(level == 3 & fold_change_ancestor < 0, max(.value_alpha2), zero)) %>%
	# 	mutate(zero = ifelse(level == 4 & fold_change_ancestor > 0, min(.value_alpha2), zero)) %>%
	# 	mutate(zero = ifelse(level == 4 & fold_change_ancestor < 0, max(.value_alpha2), zero)) %>%
	# 	
	# 	# Calculate fold chage
	# 	mutate(fold_change = .value_alpha2 - zero) %>%
	# 	mutate(significant = ((.lower_alpha2 - zero) * (.upper_alpha2 - zero)) > 0) 
		
	

#	ifelse_pipe(family == "dirichlet" | 1, ~ .x %>% get_relative_zero, ~ .x %>% mutate(zero = 0)) %>%
		
}


#------------------------------------#

run_lv_1 = function(internals,
										shards,
										level = 1,
										full_bayesian,
										approximate_posterior,
										iterations = iterations,
										sampling_iterations = sampling_iterations,
										do_regression = do_regression,
										family = family,
										.formula = .formula, model = stanmodels$ARMET_tc_fix_hierarchical){
	res1 = run_model(
		internals$reference_filtered,
		internals$mix,
		shards,
		level,
		full_bayesian,
		approximate_posterior,
		internals$prop_posterior,
		iterations = iterations,
		sampling_iterations = sampling_iterations,
		X = internals$X,
		do_regression = do_regression,
		family = family,
		cens = internals$cens,
		tree_properties = internals$tree_properties,
		Q = internals$Q,
		model = model,
		prior_survival_time = internals$prior_survival_time
	)
	
	df = res1[[1]]
	fit = res1[[2]]
	
	fit_prop_parsed = 
		fit %>%
		draws_to_tibble("prop_", "Q", "C") %>%
		filter(!grepl("_UFO|_rng", .variable))  %>%
		mutate(Q = Q %>% as.integer)
	
	draws =
		fit_prop_parsed %>%
		ungroup() %>%
		select(-.variable) %>%
		mutate(.value_relative = .value)
	
	prop =
		fit %>%
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
		left_join(df %>% distinct(Q, sample), by = "Q")	%>%
		
		# If MCMC is used check divergencies as well
		ifelse_pipe(
			!approximate_posterior,
			~ .x %>% parse_summary_check_divergence(),
			~ .x %>% parse_summary() %>% rename(.value = mean)
		) %>%
		
		# Parse
		separate(.variable, c(".variable", "level"), convert = T) %>%
		
		# Add sample information
		left_join(df %>%
								filter(`query`) %>%
								distinct(Q, sample))
	
	if (do_regression && length(parse_formula(.formula)$covariates) >0 ) 
		internals$alpha = get_alpha(fit, level, family) 

	internals$prop = prop
	internals$fit = list(fit)
	internals$df = list(df)
	internals$draws = list(draws)
	internals$prop_posterior[[1]] = fit_prop_parsed %>% group_by(.variable, Q, C) %>% prop_to_list %>% `[[` ("prop_1") 

	internals
	
	
}

