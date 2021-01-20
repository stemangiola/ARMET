clear_previous_levels = function(.data, my_level){
	# Eliminate previous results
	.data$internals$fit = .data$internals$fit[1:(my_level-1)]
	.data$internals$prop = .data$internals$prop %>% filter(level < !!my_level) 
	if(.data$internals$alpha %>% is.null %>% `!`) .data$internals$alpha = .data$internals$alpha  %>% filter(level < !!my_level)
	.data$internals$draws = .data$internals$draws[1:(my_level-1)]
	
	nodes_to_eliminate = 
		.data$proportions  %>% 
		filter(level >= !!my_level & (.variable %>% is.na %>% `!`)) %>%
		distinct(.variable) %>% 
		pull(.variable) %>%
		gsub("alpha_", "", .) %>%
		sprintf("prop_%s_prior", .)
	
	for(i in nodes_to_eliminate) {
		.data$internals$prop_posterior[[i]] = 1:ncol(.data$internals$prop_posterior[[i]]) %>% matrix(nrow=1) %>% as_tibble() %>% slice(0)
	}
	
	.data$proportions = .data$proportions  %>% filter(level < !!my_level)
	
	.data
}

#' ARMET_tc_continue
#' 
#' @description This function does inference for higher levels of the hierarchy
#' 
#' @param armet_obj An ARMEt object
#' @param level An integer
#' @param model A stan model
#' 
#' @export
ARMET_tc_continue = function(armet_obj, level, model = stanmodels$ARMET_tc_fix_hierarchical){
	
	armet_obj= clear_previous_levels(armet_obj,  level)
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
		iterations = input$iterations,
		sampling_iterations = input$sampling_iterations	,
		X = internals$X,
		do_regression = input$do_regression,
		cens = internals$cens,
		tree_properties = internals$tree_properties,
		Q = internals$Q,
		model = model,
		prior_survival_time = internals$prior_survival_time,
		sample_scaling = internals$sample_scaling,
		prior_prop = internals$prop %>% filter(level == !!level -1) %>% distinct(Q, C, .value) %>% spread(C, .value) %>% tidybulk::as_matrix(rownames = Q),
		columns_idx_including_time = internals$columns_idx_including_time,
		approximate_sampling = internals$approximate_sampling
	)
	
	df = res[[1]]
	fit = res[[2]]
	
	fit_prop_parsed = 
		fit %>%
		draws_to_tibble("prop_", "Q", "C") %>%
		filter(!grepl("_UFO|_rng|_logit", .variable))  %>%
		mutate(Q = Q %>% as.integer)
	
	draws = get_draws(fit_prop_parsed, level, internals)	
	
	prop = get_props(draws, level, df, input$approximate_posterior)	
	
	internals$prop = bind_rows(internals$prop , prop) 
	internals$fit = internals$fit %>% c(list(fit))
	internals$df = internals$df %>% c(list(df))
	internals$prop_posterior[sprintf("%s_prior", fit_prop_parsed %>% distinct(.variable) %>% pull())] = fit_prop_parsed %>% group_by(Q, C, .variable) %>% prop_to_list
	internals$draws = internals$draws %>% c(list(draws))
	
	if (input$do_regression && paste(as.character(input$.formula), collapse="")  != "~1" )
		internals$alpha = internals$alpha  %>% bind_rows( 
			get_alpha(fit, level) %>% 
				left_join(
					get_generated_quantities_standalone(fit, level, internals),
					by = c("node", "C")
					
				)
		)
	
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
				input$do_regression && paste(as.character(input$.formula), collapse="")  != "~1" ,
				~ .x %>%
					nest(proportions = -c(`Cell type category`, C, level)) %>%
					left_join(
						internals$alpha %>%	select(`Cell type category`, contains("alpha"), level, draws, rng_prop, rng_mu, .variable, one_of("Rhat")),
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
#' @importFrom tidybulk aggregate_duplicates
#' @importFrom tidybulk scale_abundance
#' @importFrom tidybulk as_matrix
#'
#' @importFrom tidyr spread
#' @importFrom tidyr gather
#' @importFrom tidyr drop_na
#' @importFrom rlang enquo
#' @importFrom tibble tibble
#' @importFrom tidybayes gather_draws
#'
#' @importFrom tidybayes gather_samples
#' @importFrom tidybayes median_qi
#'
#' @importFrom magrittr equals
#' @importFrom magrittr %$%
#'
#' @import data.tree
#'
#' @param .data A tibble
#' @param .formula A formula
#' @param .sample A column symbol
#' @param .transcript A column symbol
#' @param .abundance A column symbol
#' @param reference A tibble
#' @param approximate_posterior A boolean for variational Bayes
#' @param iterations An integer total iterations
#' @param sampling_iterations An integer. Sampling iteractions
#' @param .n_markers A tibble
#' @param do_regression A boolean
#' @param prior_survival_time An array
#' @param model A stan model
#' 
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
										iterations = 250,
										sampling_iterations = 100,
										
										.n_markers = n_markers,
										do_regression = T, 
										prior_survival_time = c(),
										model = stanmodels$ARMET_tc_fix_hierarchical,
										approximate_sampling = F) {
	 
	# At the moment is not active
	levels = 1
	cores = 4
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
	
	# Check if count is integer
	if(.data %>% select(count) %>% lapply(class) %>% unlist() %>% equals("integer") %>% `!`)
		stop(sprintf("ARMET says: the %s column must be integer as the deconvolution model is Negative Binomial", quo_name(.abundance)))
	
	# Covariate column
	if(do_regression & paste(as.character(.formula), collapse="")  != "~1"){
		formula_df = parse_formula(.formula)
		
		# Censoring column
		if(do_regression && length(formula_df$censored_column) == 1) {
			cens = .data %>% select(sample, formula_df$censored_column) %>% distinct %>% arrange(sample) %>% pull(2)
			
			# Check cens right type
			if(typeof(cens) %in% c("integer", "logical") %>% any %>% `!`) stop("ARMET says: censoring variable should be logical of integer (0,1)")
			if(length(prior_survival_time) == 0) stop("AMET says: you really need to provide third party survival time for your condition/disease")
			
			sd_survival_months = .data %>%  select(sample, formula_df$censored_value_column) %>% distinct %>% pull(formula_df$censored_value_column) %>% sd
			prior_survival_time = log1p(prior_survival_time) 
			
		}
		else{
			cens = NULL
		} 
		
		time_column = formula_df$censored_value_column 
		
		X =
			model.matrix(
				object = 	formula_df$formula_formatted,
				data = 
					.data %>% 
					select(sample, one_of(formula_df$covariates_formatted)) %>% 
					distinct %>% 
					arrange(sample) %>%
					mutate(!!as.symbol(formula_df$censored_value_column ) := log(!!as.symbol(formula_df$censored_value_column )) %>% scale())
			)
		
		columns_idx_including_time = which(grepl(time_column, colnames(X))) %>% as.array()
	}	else {
		formula_df = cens  = NULL	
		columns_idx_including_time = c()
		
		X =
			model.matrix(
				object = 	~ 1,
				data = .data %>% select(sample) %>% distinct %>% arrange(sample)
			)
	}
	
	# Do regression
	#if(length(formula_df$covariates_formatted) > 0 & (formula_df$covariates_formatted %>% is.na %>% `!`)) do_regression = T

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
					 one_of(formula_df$covariates_formatted)) %>%
		distinct() %>%
		
		# Check if data rectangular
		ifelse_pipe(
			(.) %>% check_if_data_rectangular(sample,symbol,count, type = "soft") %>% `!` &
				TRUE, #!fill_missing_values,
			~ .x %>% eliminate_sparse_transcripts(symbol)
		) %>%
		
		when(
			do_regression && length(formula_df$censored_column) == 1 ~ 
				mutate(., !!formula_df$censored_value_column := !!as.symbol(formula_df$censored_value_column) / sd_survival_months),
			~ (.)
		)
	
	mix =
		.data %>%
		select(sample, symbol, count, one_of(formula_df$covariates_formatted)) %>%
		distinct() 
	
	shards = cores #* 2
	is_level_in = shards %>% `>` (0) %>% as.integer
	
	tree = 	data.tree::Clone(ARMET::tree) 
	
	
	# Print overlap descriptive stats
	#get_overlap_descriptive_stats(mix %>% slice(1) %>% gather(symbol, count, -sample), reference)
	
	# Prepare data frames -
	# For Q query first
	# For G house keeing first
	# For GM level 1 first
	
	Q = mix %>% distinct(sample) %>% nrow
	
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
		)	
	# %>%
	# 	
	# 	# Decrease the number of house keeping used
	# 	anti_join({
	# 		mdf = (.) %>%
	# 			distinct(symbol, `house keeping`) %>%
	# 			filter(`house keeping`)
	# 		
	# 		withr::with_seed(123, 	sample_frac(mdf, 0.5)) %>%
	# 			distinct(symbol)
	# 	},
	# 	by = "symbol"
	# 	)
	# 
	tree_propeties = get_tree_properties(tree)
	
	# Find normalisation
	sample_scaling = 
		reference_filtered %>%
		filter(`house keeping`) %>% 
		select(count = lambda_log, symbol) %>%
		distinct() %>%
		mutate(count = exp(count), sample = "reference") %>% 
		tidybulk::aggregate_duplicates(sample, symbol, count, aggregation_function = median) %>%
		bind_rows(mix %>% filter(symbol %in% (reference_filtered %>% filter(`house keeping`) %>% distinct(symbol) %>% pull(symbol)))) %>%
		tidybulk::scale_abundance(sample, symbol, count, reference_sample = "reference", action ="get") %>%
		distinct(sample, multiplier)
	
	
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
			prior_survival_time = prior_survival_time,
			formula_df = formula_df,
			sample_scaling = sample_scaling,
			columns_idx_including_time = columns_idx_including_time,
			approximate_sampling = approximate_sampling
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
				do_regression && paste(as.character(.formula), collapse="")  != "~1" ,
				~ .x %>%
					nest(proportions = -c(`Cell type category`, C, level)) %>%
					left_join(
						internals$alpha %>%	select(`Cell type category`, contains("alpha"), level, draws, rng_prop, rng_mu, .variable, one_of("Rhat")),
						by = c("Cell type category", "level")
					)
			),
		
		# # Return the input itself
		input = input,
		
		# Return the fitted object
		internals = internals
	)
	
}

run_model = function(reference_filtered,
										 mix,
										 shards,
										 lv,
										 full_bayesian,
										 approximate_posterior,
										 prop_posterior,
										 iterations = 250,
										 sampling_iterations = 100,
										 X,
										 do_regression,
										 cens,
										 tree_properties,
										 Q,
										 model = stanmodels$ARMET_tc_fix_hierarchical,
										 prior_survival_time = c(),
										 sample_scaling,
										 prior_prop = matrix(1:Q)[,0],
										 columns_idx_including_time,
										 approximate_sampling) {
	
	
	# Global properties - derived by previous analyses of the whole reference dataset
	sigma_intercept = 1.3420415
	sigma_slope = -0.3386389
	sigma_sigma = 1.1720851
	lambda_mu_mu = 5.612671
	lambda_sigma = 7.131593
	
	# Non centred
	lambda_mu_prior = c(6.2, 1)
	lambda_sigma_prior =  c(3.3 , 1)
	lambda_skew_prior =  c(-2.7, 1)
	sigma_intercept_prior = c(1.9 , 0.1)
	
	# Filter on level considered
	reference_filtered = reference_filtered %>% filter(level %in% lv)
	
	df = ref_mix_format(reference_filtered, mix)
	
	G = df %>% filter(!`query`) %>% distinct(G) %>% nrow()
	GM = df %>% filter(!`house keeping`) %>% distinct(symbol) %>% nrow()
	
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
	
	# Dirichlet regression
	A = X %>% ncol
	
	# library(rstan)
	# fileConn<-file("~/.R/Makevars")
	# writeLines(c( "CXX14FLAGS += -O2","CXX14FLAGS += -DSTAN_THREADS", "CXX14FLAGS += -pthread"), fileConn)
	# close(fileConn)
	# ARMET_tc_model = rstan::stan_model("~/PhD/deconvolution/ARMET/inst/stan/ARMET_tc_fix.stan", auto_write = F)
	
	exposure_rate = 
		sample_scaling %>% 
		filter(sample %in% (y_source %>% pull(sample))) %>% 
		arrange(Q) %>% 
		pull(multiplier) 
	
	
	# Setup for exposure inference
	
	df_for_exposure = 
		df %>%
		filter(`query` &  `house keeping`)  %>%
		distinct(Q, symbol, count) %>%
		left_join(
			df %>%
				filter(!`query` & `house keeping`)  %>%
				mutate(reference_count = exp(lambda_log)) %>%
				distinct(symbol, reference_count )
		)
	
	nrow_for_exposure = nrow(df_for_exposure)
	Q_for_exposure = df_for_exposure$Q
	reference_for_exposure = df_for_exposure %>% pull(reference_count)
	counts_for_exposure = df_for_exposure %>% pull(count)
	
	init_list = list(	lambda_UFO = rep(6.2, GM)	) 
	
	ref = 
		df %>%
		
		# Eliminate the query part, not the house keeping of the query
		filter(!`query` | `house keeping`)  %>%
		
		filter(!`house keeping`) %>%
		mutate(lambda = exp(lambda_log)) %>%
		select(C, GM, lambda ) %>% 
		distinct() %>%
		arrange(C, GM) %>% 
		spread(GM, lambda) %>% 
		tidybulk::as_matrix(rownames = "C") 
	
	y = 
		y_source %>%
		select(Q, GM, count) %>% 
		distinct() %>%
		arrange(Q, GM) %>% 
		spread(GM, count) %>% 
		tidybulk::as_matrix(rownames = "Q") 
	
	max_y = max(y)
	ct_in_ancestor_level = ifelse(lv == 1, 0, tree_properties$ct_in_levels[lv-1])
	
	Sys.setenv("STAN_NUM_THREADS" = shards)
	
	if(cens %>% is.null) cens =  rep(0, Q)
	which_cens = which(cens == 1)  %>% as.array()
	which_not_cens = which(cens == 0) %>% as.array()
	how_many_cens = length(which_cens)
	
	max_unseen = ifelse(how_many_cens>0, max(X[,2]), 0 )
	if(is.null(prior_survival_time)) prior_survival_time = array(1)[0]
	spt = length(prior_survival_time)
	
	CIT = length(columns_idx_including_time)

	fit = 
		approximate_sampling %>%
		when(
			(.) ~ vb_iterative(model,
												 # rstan::stan_model("~/PhD/deconvolution/ARMET/inst/stan/ARMET_tc_fix_hierarchical.stan", auto_write = F),
												 iter = 50000,
												 tol_rel_obj = 0.0005,
												 data = prop_posterior %>% c(tree_properties),
												 init = function () init_list
			),
			
			~ 	sampling(
				model,
				#rstan::stan_model("~/PhD/deconvolution/ARMET/inst/stan/ARMET_tc_fix_hierarchical.stan", auto_write = F),
				chains = 3,
				cores = 3,
				iter = iterations,
				warmup = iterations - sampling_iterations,
				data = prop_posterior %>% c(tree_properties),
				# pars=
				# 	c("prop_1", "prop_2", "prop_3", sprintf("prop_%s", letters[1:9])) %>%
				# 	c("alpha_1", sprintf("alpha_%s", letters[1:9])) %>%
				# 	c("exposure_rate") %>%
				# 	c("lambda_UFO") %>%
				# 	c("prop_UFO") %>%
				# 	c(additional_par_to_save),
				init = function ()	init_list,
				save_warmup = FALSE
				# ,
				# control=list( adapt_delta=0.9,stepsize = 0.01,  max_treedepth =10  )
			) %>%
				{
					(.)  %>% rstan::summary() %$% summary %>% as_tibble(rownames = "par") %>% arrange(Rhat %>% desc) %>% filter(Rhat > 1.5) %>% ifelse_pipe(nrow(.) > 0, ~ .x %>% print)
					(.)
				}
		)

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

#' add_cox_test
#' 
#' @description This function adds cox regression statistics to the fit object
#' 
#' @param .data A tibble
#' @param relative A boolean
#' 
#' @export
add_cox_test = function(.data, relative = TRUE){
	
	cens_alpha = 
		.data$proportions %>% 
		select(-draws, -contains("rng")) %>%
		rename(node = .variable)  %>% 
		unnest(proportions) %>%
		censored_regression_joint(formula_df = .data$internals$formula_df, filter_how_many = Inf, relative = relative)  %>% 
		rename(.variable = node) %>%
		nest(draws_cens = -c(level, .variable  ,      C)) 
	
	.data$proportions %>%
		filter(.variable %>% is.na %>% `!`) %>%
		left_join(cens_alpha, by = c("level", "C", ".variable"))
}


#' test_differential_composition
#' 
#' @description This function performs statistical testing on the fit object
#' 
#' @param .data A tibble
#' @param credible_interval A double
#' @param cluster_CI A double
#' 
#' @export
test_differential_composition = function(.data, credible_interval = 0.90, cluster_CI = 0.55) {
       
	# x = .data$internals$formula_df$components_formatted
	# alive = .data$internals$formula_df$censored_column
	# 	
	
	cens_alpha = 
		.data$internals$formula_df$censored_formatted %>%
		when(
			length(.) == 0 ~ NULL,
			
			# Only if I have censoring
			~ .data$proportions %>% 
				select(-draws, -contains("rng")) %>%
				rename(node = .variable)  %>% 
				unnest(proportions) %>%
				censored_regression_joint(formula_df = .data$internals$formula_df, filter_how_many = Inf)  %>% 
				rename(.variable = node) %>%
				nest(draws_cens = -c(level, .variable  ,      C)) 
		)
	
	 
	.d = 
		.data$proportions %>%
		filter(.variable %>% is.na %>% `!`) %>%
		
		# If I have censoring
		when(
			!is.null(cens_alpha) ~ (.) %>% left_join(cens_alpha, by = c("level", "C", ".variable")),
			~ (.)
		) %>%
		
		cluster_posterior_slopes(credible_interval = cluster_CI) %>%
		extract_CI(credible_interval = credible_interval)
	
	dx = list()	
	
	# Level 1
	if(.d %>% filter(level ==1) %>% nrow %>% `>` (0))
		dx = 
		.d %>%
		
		filter(level ==1) %>%
		mutate(fold_change_ancestor = 0) %>%
		identify_baseline_by_clustering( ) %>%
		
		when(
			!is.null(cens_alpha) ~ 
				mutate(., significant = ((.value.lower_2 - 0) * (.value.upper_2 - 0)) > 0) %>%
				mutate(fold_change  = ifelse(significant, .value_2, 0)),
			~ mutate(., significant = ((.lower_alpha2  - 0) * (.upper_alpha2  - 0)) > 0) %>%
				mutate(fold_change  = ifelse(significant, .value_alpha2 , 0))
		) 
		
	
	# Level 2
	if(.d %>% filter(level ==2) %>% nrow %>% `>` (0))
		dx =	dx %>% bind_rows(
			.d %>%
				
				filter(level ==2) %>%
				left_join(ancestor_child, by = "Cell type category") %>%
				left_join( dx %>%	select(ancestor = `Cell type category`, fold_change_ancestor = fold_change),  by = "ancestor" )  %>%
				identify_baseline_by_clustering( ) %>%
				
				when(
					!is.null(cens_alpha) ~ 
						mutate(., significant = ((.value.lower_2 - 0) * (.value.upper_2 - 0)) > 0) %>%
						mutate(fold_change  = ifelse(significant, .value_2, 0)),
					~ mutate(., significant = ((.lower_alpha2  - 0) * (.upper_alpha2  - 0)) > 0) %>%
						mutate(fold_change  = ifelse(significant, .value_alpha2 , 0))
				) 
		) 
	
	
	# Level 3
	if(.d %>% filter(level ==3) %>% nrow %>% `>` (0))
		dx =	dx %>% bind_rows(
			.d %>%
				
				filter(level ==3) %>%
				left_join(ancestor_child, by = "Cell type category") %>%
				left_join( dx %>%	select(ancestor = `Cell type category`, fold_change_ancestor = fold_change) ,  by = "ancestor")  %>%
				identify_baseline_by_clustering( ) %>%
				
				when(
					!is.null(cens_alpha) ~ 
						mutate(., significant = ((.value.lower_2 - 0) * (.value.upper_2 - 0)) > 0) %>%
						mutate(fold_change  = ifelse(significant, .value_2, 0)),
					~ mutate(., significant = ((.lower_alpha2  - 0) * (.upper_alpha2  - 0)) > 0) %>%
						mutate(fold_change  = ifelse(significant, .value_alpha2 , 0))
				) 
		)
	
	# Level 4
	if(.d %>% filter(level ==4) %>% nrow %>% `>` (0))
		dx =	dx %>% bind_rows(
			.d %>%
				
				filter(level ==4) %>%
				left_join(ancestor_child, by = "Cell type category") %>%
				left_join( dx %>%	select(ancestor = `Cell type category`, fold_change_ancestor = fold_change) ,  by = "ancestor")  %>%
				identify_baseline_by_clustering( ) %>%
				
				when(
					!is.null(cens_alpha) ~ 
						mutate(., significant = ((.value.lower_2 - 0) * (.value.upper_2 - 0)) > 0) %>%
						mutate(fold_change  = ifelse(significant, .value_2, 0)),
					~ mutate(., significant = ((.lower_alpha2  - 0) * (.upper_alpha2  - 0)) > 0) %>%
						mutate(fold_change  = ifelse(significant, .value_alpha2 , 0))
				) 
		)
	
	dx
	

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
		cens = internals$cens,
		tree_properties = internals$tree_properties,
		Q = internals$Q,
		model = model,
		prior_survival_time = internals$prior_survival_time,
		sample_scaling = internals$sample_scaling,
		columns_idx_including_time = internals$columns_idx_including_time,
		approximate_sampling = internals$approximate_sampling
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
		draws_to_tibble("prop_1", "Q", "C") %>%
		#tidybayes::gather_draws(`prop_[1]`[Q, C], regex = T) %>%
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
		
		# If MCMC is used check divergences as well
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
	

	
	
	if (do_regression && paste(as.character(.formula), collapse="")  != "~1" ) 
		internals$alpha = 
		get_alpha(fit, level) %>% 
		left_join(
			get_generated_quantities_standalone(fit, level, internals),
			by = c("node", "C")
			
		)
	
	internals$prop = prop
	internals$fit = list(fit)
	internals$df = list(df)
	internals$draws = list(draws)
	internals$prop_posterior[[1]] = fit_prop_parsed %>% group_by(.variable, Q, C) %>% prop_to_list %>% `[[` ("prop_1") 
	
	internals
	
	
}


get_theoretical_data_disrtibution = function(fit){
	

	m2 <- rstan::stan_model(file = "inst/stan/generated_quantities_lv1.stan")
	
	
	# # If those nodes are not in fit add them otherwise generate quantities fails
	# missing_columns = 
	# 	c("1", letters[1:11]) %>%
	# 	imap(
	# 		~ sprintf("alpha_%s", .x) %>%
	# 			grep(colnames(as.matrix(fit))) %>%
	# 			when(length(.)== 0 ~ {
	# 				add_col = matrix(rep(0, nrow(as.matrix(fit)) * (tree_properties$ct_in_nodes[.y]-1) * A  ), nrow = nrow(as.matrix(fit)) )
	# 				colnames(add_col) = sprintf("alpha_%s[%s]", .x, apply(expand.grid( 1:A, 1:(tree_properties$ct_in_nodes[.y]-1)), 1, paste, collapse=","))
	# 					 
	# 					 add_col
	# 			})
	# 	) %>%
	# 	do.call(cbind,.)
	
	
}
