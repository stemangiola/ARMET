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
#' @param iterations An integer total iterations
#' @param sampling_iterations An integer. Sampling iteractions
#' @param .n_markers A tibble
#' @param do_regression A boolean
#' @param prior_survival_time An array
#' @param model A stan model
#' 
#' @rdname setup_convolved_lm
#' @name setup_convolved_lm
#'
#' @return An ARMET object
#'
#' @export
#'
setup_convolved_lm = function(.data,
										.formula = ~ 1,
										.sample = NULL,
										.transcript = NULL,
										.abundance = NULL,
										approximate_posterior = F,
										prior_survival_time = c(),
										transform_time_function = sqrt,
										...) {
	
	# At the moment is not active
	levels = 1
	full_bayesian = F
	.n_markers = n_markers
	do_regression = T
	cores = 4
	shards = cores 
	reference = NULL
	iterations = 250
	sampling_iterations = 100
	model = stanmodels$ARMET_tc_fix_hierarchical
	
	

	
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
					mutate(!!as.symbol(formula_df$censored_value_column ) := transform_time_function(!!as.symbol(formula_df$censored_value_column )) )
			)
		
		columns_idx_including_time = 
			which(grepl(time_column, colnames(X))) %>% 
			as.array() %>% 
			
			# Fix if NULL
			when(is.null(.) ~ c(), ~ (.))
		
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
		distinct(sample, multiplier) %>%
		mutate(exposure_rate = -log(multiplier)) %>%
		mutate(exposure_multiplier = exp(exposure_rate)) 
	
	
	
	# Default internals
	list(
		internals = list(
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
		approximate_posterior = approximate_posterior,
		transform_time_function = transform_time_function,
		
		shards = shards,
		levels = levels,
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

#' estimate_convoluted_lm_1
#' 
#' @description This function does inference for higher levels of the hierarchy
#' 
#' @rdname estimate_convoluted_lm
#' @name estimate_convoluted_lm_1
#' 
#' @param armet_obj An ARMET object
#' @param level An integer
#' @param model A stan model
#' 
#' @export
estimate_convoluted_lm_1 = function(armet_obj){
	
	level = 1 
	
	internals = 
		run_lv_1(
			armet_obj$internals,
			armet_obj$internals$shards,
			armet_obj$internals$levels,
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
		left_join(internals$prop) %>%
		
		# Attach alpha if regression
		ifelse_pipe(
			internals$do_regression && paste(as.character(internals$.formula), collapse="")  != "~1" ,
			~ .x %>%
				nest(proportions = -c(`Cell type category`, C, level)) %>%
				left_join(
					internals$alpha %>%	select(`Cell type category`, contains("alpha"), level, draws, rng_prop, rng_mu, .variable, one_of("Rhat")),
					by = c("Cell type category", "level")
				)
		)
	
	
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
		get_estimates(level) %>% 
		add_attr(attrib, "full_results")
}

#' estimate_convoluted_lm_2
#' 
#' @description This function does inference for higher levels of the hierarchy
#' 
#' 
#' @rdname estimate_convoluted_lm
#' @name estimate_convoluted_lm_2
#' 
#' @param armet_obj An ARMET object
#' 
#' @export
estimate_convoluted_lm_2 = function(armet_obj){
	
	level = 2
	
	attrib = attr(armet_obj, "full_results")
	
	attrib = ARMET_tc_continue(attrib, level, model = attrib$internals$model)
	
	armet_obj %>% 
		bind_rows(
			attrib$proportions %>% 
				get_estimates(level) 
		) %>% 
		
		# Add back update attributes
		add_attr(attrib, "full_results")
	
}

#' estimate_convoluted_lm_3
#' 
#' @description This function does inference for higher levels of the hierarchy
#' 
#' 
#' @rdname estimate_convoluted_lm
#' @name estimate_convoluted_lm_3
#' 
#' @param armet_obj An ARMET object
#' 
#' @export
estimate_convoluted_lm_3 = function(armet_obj){
	
	level = 3
	
	attrib = attr(armet_obj, "full_results")
	
	attrib = ARMET_tc_continue(attrib, level, model = attrib$internals$model)
	
	armet_obj %>% 
		bind_rows(
			attrib$proportions %>% 
				get_estimates(level) 
		) %>% 
		
		# Add back update attributes
		add_attr(attrib, "full_results")	
}

#' estimate_convoluted_lm_4
#' 
#' @description This function does inference for higher levels of the hierarchy
#' 
#' 
#' @rdname estimate_convoluted_lm
#' @name estimate_convoluted_lm_4
#' 
#' @param armet_obj An ARMET object
#' 
#' @export
estimate_convoluted_lm_4 = function(armet_obj){
	
	level = 4
	
	attrib = attr(armet_obj, "full_results")
	
	attrib = ARMET_tc_continue(attrib, level, model = attrib$internals$model)
	
	armet_obj %>% 
		bind_rows(
			attrib$proportions %>% 
				get_estimates(level) 
		) %>% 
		
		# Add back update attributes
		add_attr(attrib, "full_results")	
}
