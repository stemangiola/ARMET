#' ARMET-tc main
#'
#' @description This function calls the stan model.
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
#' @param prior_survival_time An array
#' @param transform_time_function transformation of the time covariate
#' @param reference A tibble
#' 
#' @rdname setup_convolved_lm
#' @name setup_convolved_lm
#'
#' @return An ARMET object
#'
#' @export
#'
#'
#'
#'

setup_convolved_lm_NON_hierarchical = function(.data,
																							 .formula = ~ 1,
																							 .sample = NULL,
																							 .transcript = NULL,
																							 .abundance = NULL,
																							 approximate_posterior = F,
																							 prior_survival_time = c(),
																							 transform_time_function = sqrt,
																							 reference = NULL,
																							 ...) {
	
	# At the moment is not active
	full_bayesian = F
	.n_markers = n_markers
	do_regression = T
	cores = 4
	shards = cores 
	iterations = 800
	sampling_iterations = 200
	model = stanmodels$ARMET_tc_fix
	
	
	
	
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
			if(length(prior_survival_time) == 0) stop("ARMET says: you really need to provide third party survival time for your condition/disease")
			
			sd_survival_months = .data %>%  select(sample, formula_df$censored_value_column) %>% distinct %>% pull(formula_df$censored_value_column) %>% sd
			prior_survival_time = transform_time_function(prior_survival_time %>% when(min(.)==0 ~ (.) + 1, (.))) 
			
			
			time_column = formula_df$censored_value_column 
			
			X =
				model.matrix(
					object = 	formula_df$formula_formatted,
					data = 
						.data %>% 
						select(sample, one_of(formula_df$covariates_formatted)) %>% 
						distinct %>% 
						arrange(sample) %>%
						mutate(!!as.symbol(formula_df$censored_value_column ) := transform_time_function(!!as.symbol(formula_df$censored_value_column )) )
				)
			
			columns_idx_including_time = 
				which(grepl(time_column, colnames(X))) %>% 
				as.array() %>% 
				
				# Fix if NULL
				when(is.null(.) ~ c(), ~ (.))
			
		}
		else{
			X =
				model.matrix(
					object = 	formula_df$formula_formatted,
					data = 
						.data %>% 
						select(sample, one_of(formula_df$covariates_formatted)) %>% 
						distinct %>% 
						arrange(sample) 
				)
			
			cens = NULL
			columns_idx_including_time = array(0)[0]
			
		} 
		
		
		
	}	else {
		formula_df = cens  = NULL	
		columns_idx_including_time = array(0)[0]
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
	
	
	
	# Print overlap descriptive stats
	#get_overlap_descriptive_stats(mix %>% slice(1) %>% gather(symbol, count, -sample), reference)
	
	# Prepare data frames -
	# For Q query first
	# For G house keeing first
	# For GM level 1 first
	
	Q = mix %>% distinct(sample) %>% nrow
	
	if(!check_if_data_rectangular(reference, cell_type, symbol, count)){
		warning("tidybulk says: the data does not have the same number of transcript per sample. The data set is not rectangular.")
		
		reference = 
			reference %>%
			
			# Filter genes common to all cell types
			add_count(symbol) %>%
			filter(n==max(n))  %>% 
			select(-n)
		
		
	}
	
	reference_filtered = reference %>% mutate(C = cell_type  %>% as.factor %>% as.numeric)
	
	# Find normalisation
	sample_scaling = 
		reference_filtered %>%
		
		mutate(sample = "reference") %>% 
		tidybulk::aggregate_duplicates(sample, symbol, count, aggregation_function = median) %>%
		bind_rows(mix) %>%
		tidybulk::identify_abundant(sample, symbol, count) %>%
		tidybulk::scale_abundance(sample, symbol, count, reference_sample = "reference", action ="get", .subset_for_scaling = .abundant) %>%
		distinct(sample, multiplier) %>%
		mutate(exposure_rate = -log(multiplier)) %>%
		mutate(exposure_multiplier = exp(exposure_rate)) 
	
	
	
	# Default internals
	list(
		internals = list(
			prop = NULL,
			fit = NULL,
			df = NULL,
			#prop_posterior = get_null_prop_posterior(tree_propeties$ct_in_nodes),
			alpha = NULL,
			Q = Q,
			reference_filtered = reference_filtered,
			mix = mix,
			X = X,
			cens = cens,
			#tree_properties = tree_propeties,
			prior_survival_time = prior_survival_time,
			formula_df = formula_df,
			sample_scaling = sample_scaling,
			columns_idx_including_time = columns_idx_including_time,
			approximate_posterior = approximate_posterior,
			transform_time_function = transform_time_function,
			
			shards = shards,
			full_bayesian = full_bayesian,
			iterations = iterations,
			sampling_iterations = sampling_iterations	,
			do_regression = do_regression,
			.formula = .formula,
			model = model
		),
		input = input
	)
	
	
}

median_qi_nest_draws_NO_hierarchical = function(d){
	# Anonymous function to add the draws to the summary
	
	left_join(
		d %>%
			#group_by(.variable,  Q,  C) %>%
			tidybayes::median_qi() %>%
			ungroup(),
		# Reattach draws as nested
		d %>%
			ungroup() %>%
			nest(.draws = c(.chain, .iteration, .draw , .value, .value_relative)),
		by = c(".variable", "Q", "sample", "C")
	)
}


# @description Parse the stan fit object and check for divergences
parse_summary_check_divergence_NO_hierarchical = function(draws) {
	draws %>%
		
		group_by(.variable,  Q, sample,  C) %>%
		
		# If not converged choose the majority chains
		mutate(converged = diptest::dip.test(`.value_relative`) %$%	`p.value` > 0.05) %>%
		
		# Anonymous function - add summary fit to converged label
	# input: tibble
	# output: tibble
	{
		left_join(
			(.) %>% select(-converged) %>% median_qi_nest_draws_NO_hierarchical(),
			(.) %>% distinct(converged),
			by = c( ".variable", "Q", "sample", "C")
		)
	} %>%
		ungroup()
}

get_generated_quantities_standalone_NO_hierarchy = function(fit, internals){
	
	
	S = internals$Q
	Q = internals$Q
	A = dim(data.frame(internals$X))[2]
	
	number_of_cell_types = internals$reference_filtered %>% distinct(cell_type) %>% nrow
	mod = stanmodels$generated_quantities
	
	fit2 = rstan::gqs(
		mod,
		draws =  as.matrix(fit)
	) 
	
	
	left_join(
		fit2 %>%
			draws_to_tibble("prop_", "Q", "C") %>%
			mutate(Q = as.integer(Q)) %>%
			mutate(.variable = gsub("_rng", "", .variable)) %>%
			separate(.variable, c("par", "node"), remove = F)  %>%
			select(-par) %>%
			nest(rng_prop = -c(node, C)) %>%
			mutate(C = 1:n()),
		
		fit2 %>%
			draws_to_tibble("mu_", "Q", "C") %>%
			mutate(Q = as.integer(Q)) %>%
			mutate(.variable = gsub("_rng", "", .variable)) %>%
			separate(.variable, c("par", "node"), remove = F)  %>%
			select(-par) %>%
			nest(rng_mu = -c(node, C)) %>%
			mutate(C = 1:n()),
		by=c("C", "node")
	)
	
	
}


#' estimate_convoluted_lm
#' 
#' @description This function does inference for higher levels of the hierarchy
#' 
#' @rdname estimate_convoluted_lm
#' @name estimate_convoluted_lm
#' 
#' @param armet_obj An ARMET object
#' 
#' @export
estimate_convoluted_lm = function(armet_obj){
	
	
	internals = 
		run_lv(
			armet_obj$internals,
			armet_obj$internals$shards,
			armet_obj$internals$full_bayesian,
			armet_obj$internals$approximate_posterior,
			iterations = armet_obj$internals$iterations,
			sampling_iterations = armet_obj$internals$sampling_iterations	,
			do_regression = armet_obj$internals$do_regression,
			.formula = armet_obj$internals$.formula,
			model = armet_obj$internals$model
		)
	

	
	proportions =	
		armet_obj$input$.data %>%
		select(c(sample, (.) %>% get_specific_annotation_columns(sample))) %>%
		distinct() %>%
		left_join(internals$prop, by="sample") %>%
		
		# Attach alpha if regression
		ifelse_pipe(
			internals$do_regression, # && paste(as.character(internals$.formula), collapse="")  != "~1" ,
			~ .x %>%
				nest(proportions = -c( C)) %>%
				left_join(
					internals$alpha %>%	select( C, contains("alpha"), draws, rng_prop, rng_mu, .variable, one_of("Rhat")),
					by = c("C")
				)
		) %>% 
		
		# Add cell type
		left_join(internals$reference_filtered %>% distinct(cell_type, C), by="C") %>% 
		select(-C) %>% 
		select(cell_type, everything())
	
	attrib = 
		list(
			# Matrix of proportions
			proportions = proportions,
			
			# # Return the input itself
			input = armet_obj$input,
			
			# Return the fitted object
			internals = internals
		)
	
	proportions %>% 
		get_estimates_NO_hierarchy(X = attrib$internals$X) %>% 
		add_attr(attrib, "full_results")
}

run_lv = function(internals,
									shards,
									full_bayesian,
									approximate_posterior,
									iterations = iterations,
									sampling_iterations = sampling_iterations,
									do_regression = do_regression,
									.formula = .formula, model = stanmodels$ARMET_tc_fix){
	
	res1 = run_model_no_hierarchy( 
		internals$reference_filtered,
		internals$mix,
		shards,
		full_bayesian,
		approximate_posterior,
		internals$prop_posterior,
		iterations = iterations,
		sampling_iterations = sampling_iterations,
		X = internals$X,
		do_regression = do_regression,
		cens = internals$cens,
		Q = internals$Q,
		model = model,
		prior_survival_time = internals$prior_survival_time,
		sample_scaling = internals$sample_scaling,
		columns_idx_including_time = internals$columns_idx_including_time
		
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
		
		# add sample annotation
		left_join(df %>% distinct(Q, sample), by = "Q")	%>%
		
		# If MCMC is used check divergences as well
		ifelse_pipe(
			!approximate_posterior,
			~ .x %>% parse_summary_check_divergence_NO_hierarchical(),
			~ .x %>% parse_summary() %>% rename(.value = mean)
		) %>%
		
		# Parse
		separate(.variable, c(".variable", "level"), convert = T) %>%
		
		# Add sample information
		left_join(df %>%
								filter(`query`) %>%
								distinct(Q, sample))
	
	
	
	
	if (do_regression) # && paste(as.character(.formula), collapse="")  != "~1" ) 
		internals$alpha = 
		get_alpha_NO_hierarchy(fit) %>% 
		left_join(
			get_generated_quantities_standalone_NO_hierarchy(fit, internals),
			by = c("node", "C")
			
		)
	
	
	internals$prop = prop
	internals$fit = list(fit)
	internals$df = list(df)
	internals$draws = list(draws)
	internals$prop_posterior[[1]] = fit_prop_parsed %>% group_by(.variable, Q, C) %>% prop_to_list %>% `[[` ("prop_1") 
	
	internals
	
	
}

get_estimates_NO_hierarchy = function(.data, X) {
	
	.data %>% 
		filter(.variable %>% is.na %>% `!`) %>%
		select(cell_type, draws) %>% 
		mutate(regression = map(draws,
														~ .x %>%
															group_by(A) %>%
															summarise(.median = median(.value), .sd = sd(.value)) %>% 
															#tidybayes::median_qi(.width = credible_interval) %>%
															
															left_join(tibble(A=1:ncol(X), A_name = colnames(X)) ,  by = "A") %>% 
															select(-A) %>% 
															
															pivot_wider(
																names_from = A_name,
																values_from = c(.median, .sd)
															))) %>% 
		select(-draws) %>% 
		unnest(regression)
	
}

run_model_no_hierarchy = function(reference_filtered,
																	mix,
																	shards,
																	full_bayesian,
																	approximate_posterior,
																	prop_posterior,
																	iterations = 250,
																	sampling_iterations = 100,
																	X,
																	do_regression,
																	cens,
																	Q,
																	model = stanmodels$ARMET_tc_fix,
																	prior_survival_time = c(),
																	sample_scaling,
																	prior_prop = matrix(1:Q)[,0, drop=FALSE],
																	columns_idx_including_time) {
	
	
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
	my_genes = reference_filtered %>% filter(is_marker) %>% pull(symbol) %>% unique()
	reference_filtered = reference_filtered %>% filter(symbol %in% my_genes)
	
	# Number of cell types
	number_of_cell_types = reference_filtered %>% distinct(cell_type) %>% nrow()
	
	# Format
	df = ref_mix_format(reference_filtered, mix)
	
	GM = df  %>% distinct(symbol) %>% nrow()
	
	y_source =
		df %>%
		filter(`query`) %>%
		select(S, Q, symbol, count, GM, sample) 
	
	# Dirichlet regression
	A = X %>% ncol
	
	# library(rstan)
	# fileConn<-file("~/.R/Makevars")
	# writeLines(c( "CXX14FLAGS += -O2","CXX14FLAGS += -DSTAN_THREADS", "CXX14FLAGS += -pthread"), fileConn)
	# close(fileConn)
	# ARMET_tc_model = rstan::stan_model("~/PhD/deconvolution/ARMET/inst/stan/ARMET_tc_fix.stan", auto_write = F)
	
	exposure_multiplier = 
		sample_scaling %>% 
		filter(sample %in% (y_source %>% pull(sample))) %>% 
		arrange(sample) %>% 
		pull(exposure_multiplier) %>%
		as.array()

	
	init_list = list(	lambda_UFO = rep(6.2, GM)	) 
	
	ref = 
		df %>%
		
		# Eliminate the query part, not the house keeping of the query
		filter(!`query`)  %>%
		
		select(C, GM, count ) %>% 
		distinct() %>%
		arrange(C, GM) %>% 
		spread(GM, count) %>% 
		tidybulk::as_matrix(rownames = "C") 
	
	y = 
		y_source %>%
		select(Q, GM, count) %>% 
		distinct() %>%
		arrange(Q, GM) %>% 
		spread(GM, count) %>% 
		tidybulk::as_matrix(rownames = "Q") 
	
	max_y = max(y)

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
		approximate_posterior %>%
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
				#data = prop_posterior %>% c(tree_properties),
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
	
	list(df, fit)
	
}