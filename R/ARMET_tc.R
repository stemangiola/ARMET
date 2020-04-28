#' ARMET_tc_continue
#' 
#' @description This function
#' 
#' @export
ARMET_tc_continue = function(armet_obj, levels, model = stanmodels$ARMET_tc_fix_hierarchical){
	
	run_lv = switch(	levels,run_lv_1,	run_lv_2,	run_lv_3,	run_lv_4)
	
	
	internals = 
		run_lv(
			armet_obj$internals,
			shards = armet_obj$input$cores,
			level = levels,
			full_bayesian = armet_obj$input$full_bayesian,
			approximate_posterior = armet_obj$input$approximate_posterior,
			iterations = armet_obj$input$iterations,
			sampling_iterations = armet_obj$input$sampling_iterations	,
			do_regression = armet_obj$input$do_regression,
			family = armet_obj$input$family,
			.formula = armet_obj$input$.formula,
			model = model
		)
	
	# Return
	list(
		# Matrix of proportions
		proportions =	
			armet_obj$input$.data %>%
			select(c(sample, (.) %>% get_specific_annotation_columns(sample))) %>%
			distinct() %>%
			left_join(internals$prop) %>%
			
			# Attach alpha if regression
			ifelse_pipe(
				armet_obj$input$do_regression && length(parse_formula(armet_obj$input$.formula)$covariates) >0 ,
				~ .x %>%
					nest(proportions = -c(`Cell type category`, C, level)) %>%
					left_join(
						internals$alpha %>%	select(`Cell type category`, contains("alpha"), zero, level, draws, rng, .variable, Rhat),
						by = c("Cell type category", "level")
					)
			),
		
		# # Return the input itself
		input = armet_obj$input,
		
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
										levels = 3,
										.n_markers = n_markers ,
										do_regression = T, 
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
	
	# Checkif count is integer
	if(.data %>% select(count) %>% lapply(class) %>% unlist() %>% equals("integer") %>% `!`)
		stop(sprintf("ARMET says: the %s column must be integer as the deconvolution model is Negative Binomial", quo_name(.abundance)))

	# Check family
	if(family %in% c("dirichlet", "beta") %>% any %>% `!`) stop("ARMET says: Please choose between dirichlet or beta families")

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
		
		sd_survival_months = 29.3
		
		df_for_edgeR = df_for_edgeR %>% mutate(!!.cens_value_column := !!as.symbol(.cens_value_column) / sd_survival_months)
		#surv_prior = .data %>% filter(!(!!as.symbol(.cens_column))) %>% select(sample, .cens_column, cov_columns) %>% distinct() %>% pull(cov_columns[1]) %>% gamma_alpha_beta ()

		# Check cens right type
		if(typeof(cens) %in% c("integer", "logical") %>% any %>% `!`) stop("ARMET says: censoring variable should be logical of integer (0,1)")
	}
	else{
		cens = NULL
		#surv_prior = c()
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
		select(sample,
					 symbol,
					 count,
					 one_of(parse_formula(.formula)$covariates)) %>%
		distinct() %>%
		spread(symbol, count)

	shards = cores #* 2
	is_level_in = shards %>% `>` (0) %>% as.integer

	tree = 	data.tree::Clone(ARMET::tree) 
	# %>%	{
	# 	# Filter selected levels
	# 	data.tree::Prune(., function(x)
	# 		x$level <= max(levels_in_the_tree) + 1)
	# 
	# 	# Filter if not in referenc
	# 	#data.tree::Prune(., function(x) ( x$name %in% (ARMET::ARMET_ref %>% distinct(`Cell type category`) %>% pull(1) %>% as.character) ))
	# 	.
	# }


	
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

	run_lv = switch(	levels,run_lv_1,	run_lv_2,	run_lv_3,	run_lv_4)
	
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
			tree_properties = tree_propeties
		) 
	
	internals = 
		run_lv(
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
						internals$alpha %>%	select(`Cell type category`, contains("alpha"), zero, level, draws, rng, .variable, Rhat),
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
										 model = stanmodels$ARMET_tc_fix_hierarchical) {
	
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
	which_cens = which(cens == 1)
	which_not_cens = which(cens == 0)
	how_many_cens = length(which_cens)

	max_unseen = ifelse(how_many_cens>0, max(X[,2]), 0 )


	
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

# lower_triangular = function(.data){
# 	
# 	levs = .data$`Cell type category_1` %>% levels
# 	
# 	.data %>%
# 		select(`Cell type category_1`, `Cell type category_2`,    prob) %>%
# 		spread(`Cell type category_2` ,   prob) %>% 
# 		as_matrix(rownames = "Cell type category_1") %>%
# 		
# 		# Drop upper triangular
# 		{ ma = (.); ma[lower.tri(ma)] <- NA; ma} %>% 
# 		
# 		as_tibble(rownames = "Cell type category_1") %>% 
# 		gather(`Cell type category_2`, prob, -`Cell type category_1`) %>% 
# 		mutate(
# 			`Cell type category_1` = factor(`Cell type category_1`, levels = levs), 
# 			`Cell type category_2` = factor(`Cell type category_2`, levels = levs), 
# 		) %>%
# 		drop_na
# }

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
				permute_nest(
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
						mutate(frac_up = map2_dbl(data, med, ~ (.x$.value  > .y) %>% sum %>% divide_by(nrow(.x)))) %>% 
						mutate(frac_down = map2_dbl(data, med, ~ (.x$.value  < .y) %>% sum %>% divide_by(nrow(.x)))) %>%
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

plot_heatmap = function(.data){
	
	.data %>%
		get_signatures %>%
		filter(A==2, .variable %in% c("alpha_1", "alpha_a")) %>%
		lower_triangular %>%
		separate(`Cell type category_1`, c("l1"), sep="_", extra = "drop", remove = F) %>%
		separate(`Cell type category_2`, c("l2"), sep="_", extra = "drop", remove = F) %>%
		
		mutate(
			label = paste(
				substr(`l2`, start = 0, stop = 3),
				substr(`l1`, start = 0, stop = 3),
				sep = "\n"
			)
		) %>%
		ggplot(aes( `Cell type category_1`, `Cell type category_2`, fill=prob, label=label)) +
		geom_tile() +
		geom_text(angle = 45) +
		#facet_wrap(~.variable, scale="free") +
		scale_fill_distiller(palette = "Spectral", na.value="transparent", limits = c(-1,1) ) +
		coord_fixed() +
		scale_x_discrete(position = "top") +
		theme(
			plot.background = element_rect(fill = "transparent",colour = NA),
			panel.grid=element_blank(),
			panel.background=element_blank(),
			panel.border = element_blank(),
			plot.margin = unit(c(0, 0, 0, 0), "npc"),
			axis.title.x = element_blank(),
			axis.title.y = element_blank(),
			axis.text = element_blank()
			#,
			#axis.text.x = element_text(angle = 90, hjust = 0)
		)
	
	# plot_2 <- ARMET_TCGA_result_hierarchical_beta %>%
	# 	get_signatures %>%
	# 	filter(A==2) %>%
	# 	nest(data = -.variable) %>%
	# 	mutate(lower_tri = map(data, ~.x %>% lower_triangular)) %>%
	# 	mutate(p = map(
	# 		lower_tri,
	# 		~.x %>%
	# 			ggplot(aes( `Cell type category_1`, `Cell type category_2`, fill=prob)) +
	# 			geom_tile() +
	# 			#facet_wrap(~.variable, scale="free") +
	# 			scale_fill_distiller(palette = "Spectral", na.value="transparent", limits = c(-1,1) ) +
	# 			coord_fixed() +
	# 			scale_x_discrete(position = "top") +
	# 			theme(
	# 				plot.background = element_rect(fill = "transparent",colour = NA),
	# 				panel.grid=element_blank(),
	# 				panel.background=element_blank(),
	# 				panel.border = element_blank(),
	# 				plot.margin = unit(c(0, 0, 0, 0), "npc"),
	# 				axis.title.x = element_blank(),
	# 				axis.title.y = element_blank(),
	# 				axis.text.x = element_text(angle = 90, hjust = 0)
	# 			)
	# 	))
	
}



#' @export
test_differential_composition = function(.data, credible_interval = 0.90, cluster_CI = 0.55) {


	.d = 
		.data$proportions %>%
		filter(.variable %>% is.na %>% `!`) %>%
		select(-one_of("zero")) %>%
		cluster_posterior_slopes(credible_interval = cluster_CI) %>%
		extract_CI
	
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
		left_join(ARMET::tree %>% get_ancesotr_child) %>%
		left_join( dx %>%	select(ancestor = `Cell type category`, fold_change_ancestor = fold_change) )  %>%
		identify_baseline_by_clustering( ) %>%
		
		mutate(significant = ((.lower_alpha2 - zero) * (.upper_alpha2 - zero)) > 0) %>%
		mutate(fold_change  = ifelse(significant, .value_alpha2, 0))
		) 
	
	
	# Level 3
	if(.d %>% filter(level ==3) %>% nrow %>% `>` (0))
		dx =	dx %>% bind_rows(
		.d %>%
		
		filter(level ==3) %>%
		left_join(ARMET::tree %>% get_ancesotr_child) %>%
		left_join( dx %>%	select(ancestor = `Cell type category`, fold_change_ancestor = fold_change) )  %>%
		identify_baseline_by_clustering( ) %>%
		
		mutate(significant = ((.lower_alpha2 - zero) * (.upper_alpha2 - zero)) > 0) %>%
		mutate(fold_change  = ifelse(significant, .value_alpha2, 0))
		)
	
	# Level 4
	if(.d %>% filter(level ==4) %>% nrow %>% `>` (0))
		dx =	dx %>% bind_rows(
			.d %>%
				
				filter(level ==4) %>%
				left_join(ARMET::tree %>% get_ancesotr_child) %>%
				left_join( dx %>%	select(ancestor = `Cell type category`, fold_change_ancestor = fold_change) )  %>%
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
	# 	left_join(ARMET::tree %>% get_ancesotr_child)	%>%
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

#' Create polar plot of results
#' @rdname plot_polar
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

	xx  = 
	.data %>%
		calculate_x_for_polar %>%
		
		# Add angle text
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
		unnest(summarised_proportion) %>%
		#mutate(median_proportion = median_proportion / max(median_proportion)) %>%
		
		distinct() %>%
		left_join(
			ct_names_polar, by=c("Cell type category" = "name")
		)

	my_rescale = function(x) { as.character( format(round( x * xx %>% pull(median_proportion) %>% max, 2), nsmall = 2)) }
	cusotm_root_trans = function() scales::trans_new("cusotm_root",function(x) x^(1/4), function(x) x^(4))
	DTC_scale =  xx %>% pull(Estimate) %>% abs() %>% { switch( (!all(is.na(.))) + 1, 1, .) } %>% max(na.rm = T)
	y_text_out = 8.5
	y_text_in = 1.4641

	# Size outer circles
	soc = (c(0, 0.2, 0.4, 0.6, 0.8)*0.5 + 1)^4  

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
		annotate("rect", xmin=0, xmax=1, ymin=soc[1], ymax= soc[2], fill="grey95") +
		annotate("rect", xmin=0, xmax=1, ymin=soc[2], ymax=soc[3], fill="grey90") +
		annotate("rect", xmin=0, xmax=1, ymin=soc[3], ymax=soc[4], fill="grey87") +
		annotate("rect", xmin=0, xmax=1, ymin=soc[4], ymax=soc[5], fill="grey83") +
		
		# Lv 1
		geom_bar(
			data = xx %>% filter(level == 1),
			aes(width = leafCount_norm, y = `median_proportion`), # / leafCount_norm),
			color = "grey20",  stat = "identity"
		) +
		geom_errorbar(
			data = xx %>% filter(level == 1), # %>% mutate(x = ifelse(`Cell type category`=="immune_cell", 0.5, x)),
			aes(width = 0 , ymin=`median_proportion`, ymax=mean(soc[1:2])), # / leafCount_norm),
			color = "grey20",  stat = "identity"
		) +
		geom_text(
			data =
				xx %>%
				filter(level == 1) %>%
				# mutate(x = ifelse(`Cell type category`=="immune_cell", 0.5, x)) %>%
				# mutate(angle = ifelse(`Cell type category`=="immune_cell", -0, angle)) %>%
				mutate(`formatted` = sub("\\s+$", "", `formatted`)),
			aes(label=`formatted`, y = mean(soc[1:2]), angle= angle  ) ,size =size_geom_text ) +
		#scale_x_continuous(labels = xx %>% filter(level == 1) %>% pull(`formatted`), breaks = xx %>% filter(level == 1) %>% pull(leafCount_norm_cum) - 0.5 * xx %>% filter(level == 1) %>% pull(leafCount_norm)) +
		
		# Lv 2
		geom_bar(
			data = xx %>% filter(level == 2),
			aes(width = leafCount_norm, y = `median_proportion` ), # / leafCount_norm),
			color = "grey20",  stat = "identity"
		) +
		geom_errorbar(
			data = xx %>% filter(level == 2),
			aes(width = 0 , ymin=`median_proportion`, ymax=mean(soc[2:3])), # / leafCount_norm),
			color = "grey20",  stat = "identity",linetype="dotted"
		) +
		geom_text(
			data =
				xx %>%
				filter(level == 2) %>%
				mutate(`formatted` = sub("\\s+$", "", `formatted`)),
			aes(label=`formatted`, y = mean(soc[2:3]), angle= angle) ,size =size_geom_text) +
		#scale_x_continuous(labels = xx %>% filter(level == 2) %>% pull(`formatted`), breaks = xx %>% filter(level == 2) %>% pull(leafCount_norm_cum) - 0.5 * xx %>% filter(level == 2) %>% pull(leafCount_norm)) +
		
		# Lv 3
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
					aes(width = 0 , ymin=`median_proportion`, ymax=mean(soc[3:4])), # / leafCount_norm),
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
						mutate(`formatted` = sub("\\s+$", "", `formatted`)),
					aes(label=`formatted`, y = mean(soc[3:4]) , angle= angle) ,size =size_geom_text
				),

				# If there are NOT cell types
				geom_errorbar(ymin=0, ymax=0)
			)
		} +
		
		
		# Lv 4
		geom_bar(
			data = xx %>% filter(level == 4),
			aes(width = leafCount_norm, y = `median_proportion` ), # / leafCount_norm),
			color = "grey20", stat = "identity"
		)  +
		{
			# Make plotting robust if no level 4 cell types were detected
			switch(
				(!  xx %>% filter(level == 4) %>% filter(`median_proportion`>prop_filter | !is.na(Estimate)) %>% nrow() > 0) + 1,
				
				# If there are cell types
				geom_errorbar(
					data = xx %>% filter(level == 4) %>% filter(`median_proportion`>prop_filter | !is.na(Estimate)),
					aes(width = 0 , ymin=`median_proportion`, ymax=mean(soc[4:5])), # / leafCount_norm),
					color = "grey20",  stat = "identity",linetype="dotted"
				) ,
				
				# If there are NOT cell types
				geom_errorbar(ymin=0, ymax=0)
			)
		} +
		{
			# Make plotting robust if no level 4 cell types were detected
			switch(
				(!  xx %>% filter(level == 4) %>% filter(`median_proportion`>prop_filter | !is.na(Estimate)) %>% nrow() > 0) + 1,
				
				# If there are cell types
				
				geom_text(
					data =
						xx %>%
						filter(level == 4)  %>%
						filter(
							`median_proportion`>prop_filter |
								!is.na(Estimate)
						) %>%
						mutate(`formatted` = sub("\\s+$", "", `formatted`)),
					aes(label=`formatted`, y = mean(soc[4:5]) , angle= angle) ,size =size_geom_text
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

#' 
#' #' Create polar plot of results
#' @rdname plot_scatter
#'
#' Prints a report of the hipothesis testing
#'
#' @import ggplot2
#' @import tibble
#' @import dplyr
#'
#'
#' @param ARMET-tc object
#'
#' @return a ggplot
#'
#' @export
#' 
plot_scatter = function(.data){
	
	data_CI = 
		.data$proportions %>%
		select(-proportions, -.variable) %>%
		#unnest(draws) %>%
		select( -draws) %>%
		unnest(rng)  %>% group_by(zero, C, Q, .variable, level, `Cell type category`, Rhat) %>% tidybayes::median_qi() %>% ungroup() %>%
		left_join(
			.data$proportions %>%
				select(-draws, -.variable) %>%
				unnest(proportions) %>% distinct(Q,  DFS_MONTHS)
		) %>%
		select(level,  `Cell type category`, Q, .upper, .lower)
	

		inferred_x = 
			map2_dfr(
			.data$internals$fit,
			1:length(.data$internals$fit),
			~ .x %>%
				tidybayes::gather_draws(X_[Q, A]) %>% 
				tidybayes::median_qi() %>% 
				ungroup() %>% 
				filter(A==2) %>% 
				mutate(level = .y)
		) %>% rename(inferred_x = .value)
	
	.data$proportions %>%
		select(-draws, -.variable) %>%
		unnest(proportions) %>%
		left_join(data_CI) %>%
		left_join(inferred_x %>% select(level, Q, inferred_x)) %>%
		ggplot(aes(x = inferred_x, y = `.value_relative`, label=sample)) +
		# geom_point(
		# 	data =
		# 		res$proportions %>% filter(level ==1) %>%
		# 		select(-proportions, -.variable) %>%
		# 		#unnest(draws) %>%
		# 		select( -draws) %>%
		# 		unnest(rng) %>%
		# left_join(
		# 	res$proportions %>%
		# 		filter(level ==1) %>%
		# 		select(-draws, -.variable) %>%
	# 		unnest(proportions) %>% distinct(Q,  DFS_MONTHS)
	# ),
	# 	aes(x = DFS_MONTHS, y = .value),
	# 	color="black",
	# 	alpha=0.1) +
	geom_errorbar( aes(ymin=.lower, ymax=.upper) ) +
	geom_point(aes(color = alive)) +
	facet_wrap(~`Cell type category`, scale="free")
}

#' This is a generalisation of ifelse that acceots an object and return an objects
#'
#' @import ggplot2
#'
#' @export
plot_markers  = function(result, level, S = NULL, cores = 20){
	
	library(multidplyr)
	
	result$signatures[[level]] %>%
		filter(!query & !`house keeping`) %>%
		distinct(`Cell type category`, C, level, G, GM, symbol, lambda_log, sigma_inv_log) %>%
		{ print(1); Sys.time(); (.) } %>%
		
		# Add divergence
		left_join(
			result$fit[[level]] %>% rstan::summary() %$% summary %>% as_tibble(rownames = "par") %>% filter(Rhat > 1.6) %>%
				filter(grepl("^lambda_log", par)) %>%
				separate(par, c("par", "G"), sep="\\[|\\]", extra = "drop") %>%
				distinct(G) %>% mutate(G = G %>% as.integer) %>%
				mutate(converged = F)
		) %>%
		mutate(converged = ifelse(converged %>% is.na, T, converged)) %>%
		{ print(2); Sys.time(); (.) } %>%
		
		# If inferred replace lambda_log and sigma_inv_log
		ifelse_pipe(
			result$input$full_bayesian,
			~ .x %>% select(-lambda_log, -sigma_inv_log) %>%
				left_join(
					result$fit[[level]] %>% tidybayes::spread_draws(lambda_log[G], sigma_inv_log[G]) %>% ungroup() %>%
						rename(lambda_log = lambda_log,sigma_inv_log = sigma_inv_log)
				)
		) %>%
		{ print(3); Sys.time(); (.) } %>%
		
		left_join(
			ARMET::ARMET_ref %>% distinct(symbol, ct1, ct2, level)
		) %>%
		
		# Add proportions
		left_join(
			result$proportions %>% select(level, Q, sample, .draws, `Cell type category`) %>% unnest(.draws)
		) %>%
		
		# add expsure
		left_join(
			result$fit[[level]] %>% tidybayes::spread_draws(exposure_rate[S]) %>% ungroup() %>% rename(Q = S)
		) %>%
		
		# Filter by sample
		ifelse_pipe(
			S %>% is.null %>% `!`,
			~ .x %>% filter(Q == S)
		)	 %>%
		
		# Calculate sum
		mutate(
			lambda_exp = lambda_log %>% exp,
			sigma_exp = 1 / exp(sigma_inv_log)
		) %>%
		
		{ print(4); Sys.time(); (.) } %>%
		
		# Filter just first 30 draws
		inner_join( (.) %>% distinct(.draw) %>% sample_n(30) ) %>%
		
		do_parallel_start(cores, "symbol") %>%
		do({
			
			`%>%` = magrittr::`%>%`
			
			sum_NB = function(lambda, sigma, prop){
				
				prop_mat = matrix(prop, nrow=1)
				lambda_mat = matrix(lambda, ncol = 1)
				sigma_mat = matrix(sigma, ncol = 1)
				
				lambda_sum = prop_mat %*% lambda_mat;
				sigma_sum =
					lambda_sum^2 /
					(
						prop_mat %*%
							(
								lambda_mat^2 /
									sigma_mat
							)
					) ;
				
				
				c(lambda_sum, sigma_sum)
				
			}
			
			(.) %>%
				tidyr::nest(data_for_sum = -c(level, symbol, GM, converged, ct1    ,     ct2   ,     sample   ,   exposure_rate,   .chain, .iteration ,.draw )) %>%
				dplyr::mutate(.sum = purrr::map(
					data_for_sum,
					~
						sum_NB(.x$lambda_exp, .x$sigma_exp, .x$.value) %>%
						as.matrix %>%
						t %>%
						tibble::as_tibble() %>%
						setNames(c("lambda_sum", "sigma_sum"))
					
				))
			
		}) %>%
		do_parallel_end() %>%
		
		{ print(5); Sys.time(); (.) } %>%
		
		select(-data_for_sum) %>%
		unnest(.sum) %>%
		{ print(6); Sys.time(); (.) } %>%
		
		# normalise
		mutate(lambda_sum = lambda_sum * exp(exposure_rate)) %>%
		
		# Calculate generated quantities
		mutate(counts_inferred = rnbinom(n(), mu = lambda_sum, size = sigma_sum)) %>%
		{ print(7); Sys.time(); (.) } %>%
		
		# Summarise
		group_by(level, symbol, GM, converged, ct1, ct2, sample, .chain) %>%
		summarize(.mean = mean(counts_inferred),
							.sd = sd(counts_inferred),
							.q025 = quantile(counts_inferred, probs = .025),
							.q25 = quantile(counts_inferred, probs = .25),
							.q50 = quantile(counts_inferred, probs = .5),
							.q75 = quantile(counts_inferred, probs = .75),
							.q97.5 = quantile(counts_inferred, probs = .975)
		) %>%
		{ print(8); Sys.time(); (.) } %>%
		
		# Add counts
		left_join(	result$input$mix %>%	gather(symbol, count, -sample) ) %>%
		{ print(9); Sys.time(); (.) } %>%
		
		{ ((.) %>%	ggplot(aes(x = count+1, y=.q50+1, color=ct1, shape = converged, GM = GM)) +
			 	geom_abline() +
			 	geom_errorbar(aes(ymin = .q025 + 1, ymax=.q97.5 + 1), alpha=0.5) +
			 	geom_point() +
			 	scale_y_log10() +
			 	scale_x_log10() +
			 	facet_grid(converged ~.chain) +
			 	my_theme) %>% print
			
			(.)
		}
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
		model = model
	)
	
	df1 = res1[[1]]
	fit1 = res1[[2]]
	
	fit_prop_parsed = fit1 %>%
		tidybayes::gather_draws(prop_1[Q, C])
	
	draws_1 =
		fit_prop_parsed %>%
		ungroup() %>%
		select(-.variable) %>%
		mutate(.value_relative = .value)
	
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
	
	if (do_regression && length(parse_formula(.formula)$covariates) >0 ) {
		alpha_1 =
			fit1 %>%
			tidybayes::gather_draws(`alpha_[1]`[A, C], regex = T) %>%
			ungroup() %>%
			
			# rebuild the last component sum-to-zero
			ifelse_pipe(family == "dirichlet" | 1, ~ .x %>% rebuild_last_component_sum_to_zero) %>%
			
			# Calculate relative 0 because of dirichlet relativity
			ifelse_pipe(family == "dirichlet" | 1, ~ .x %>% get_relative_zero, ~ .x %>% mutate(zero = 0)) %>%
			
			arrange(.chain, .iteration, .draw,     A) %>%
			
			nest(draws = -c(C, .variable, zero)) %>%
			
			# Attach convergence information
			left_join(
				rstan::summary(fit1, pars=c("alpha_1"))$summary %>% 
					as_tibble(rownames=".variable") %>% 
					separate(.variable, c(".variable", "A", "C"), sep="\\[|\\]|,", extra = "drop") %>% 
					mutate(A = as.integer(A), C = as.integer(C)) %>% filter(A == 2) %>% 
					select(.variable, C, Rhat),
				by = c(".variable", "C")
			) %>%
			
			left_join(
				tree %>%
					ToDataFrameTree("name", "level", sprintf("C%s", 1)) %>%
					filter(level == 1 + 1) %>%
					arrange(C1) %>%
					mutate(C = 1:n()) %>%
					select(name, C)  %>%
					rename(`Cell type category` = name)
			)
		
		alpha_1 = alpha_1 %>%
			separate(.variable, c("par", "node"), remove = F) %>%
			left_join(
				fit1 %>% tidybayes::gather_draws(`prop_[1a-z]_rng`[Q, C], regex = T) %>% ungroup() %>% mutate(.variable = gsub("_rng", "", .variable)) %>% separate(.variable, c("par", "node"), remove = F)  %>% select(-par) %>% nest(rng = -c(node, C)) 
			)
		
		internals$alpha = alpha_1 %>% mutate(level = 1)
		
		
	}
	
	internals$prop = prop_1
	internals$fit = list(fit1)
	internals$df = list(df1)
	internals$draws = list(draws_1)
	internals$prop_posterior[[1]] = fit_prop_parsed %>% prop_to_list %>% `[[` ("prop_1") 

	internals
	
	
}




#------------------------------------#
run_lv_2 = function(internals,
										shards,
										level = 2,
										full_bayesian,
										approximate_posterior,
										iterations = iterations,
										sampling_iterations = sampling_iterations,
										do_regression = do_regression,
										family = family,
										.formula = .formula, model = stanmodels$ARMET_tc_fix_hierarchical){
	res2 = run_model(
		internals$reference_filtered,
		internals$mix,
		shards,
		level,
		full_bayesian,
		approximate_posterior,
		internals$prop_posterior,
		draws_to_exposure(internals$fit[[1]])	,
		iterations = iterations,
		sampling_iterations = sampling_iterations,
		X = internals$X,
		do_regression = do_regression,
		family = family,
		cens = internals$cens,
		tree_properties = internals$tree_properties,
		Q = internals$Q,
		model = model
	)
	
	
	
	df2 = res2[[1]]
	fit2 = res2[[2]]
	
	fit_prop_parsed = fit2 %>%
		tidybayes::gather_draws(prop_a[Q, C])
		
	draws_a =
		fit_prop_parsed %>%
		ungroup() %>%
		select(-.variable)
	
	draws_2 =
		internals$draws[[1]] %>%
		rename(C1 = C) %>%
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
						C1 = internals$tree_properties$parents_lv2
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
	
	if (do_regression && length(parse_formula(.formula)$covariates) >0 ){
		alpha_2 =
			fit2 %>%
			tidybayes::gather_draws(`alpha_[a]`[A, C], regex = T) %>%
			ungroup() %>%
			
			# rebuild the last component sum-to-zero
			ifelse_pipe(family == "dirichlet" | 1, ~ .x %>% rebuild_last_component_sum_to_zero) %>%
			
			# Calculate relative 0 because of dirichlet relativity
			ifelse_pipe(family == "dirichlet" | 1, ~ .x %>% get_relative_zero, ~ .x %>% mutate(zero = 0)) %>%
			
			arrange(.chain, .iteration, .draw,     A) %>%
			
			nest(draws = -c(C, .variable, zero)) %>%
			
			# Attach convergence information
			left_join(
				rstan::summary(fit2, pars=c("alpha_a"))$summary %>% 
					as_tibble(rownames=".variable") %>% 
					separate(.variable, c(".variable", "A", "C"), sep="\\[|\\]|,", extra = "drop") %>% 
					mutate(A = as.integer(A), C = as.integer(C)) %>% filter(A == 2) %>% 
					select(.variable, C, Rhat),
				by = c(".variable", "C")
			) %>%
			
			left_join(
				# FOR HIERARHICAL
				tree %>%
					ToDataFrameTree("name", "level", sprintf("C%s", 2)) %>%
					filter(level == 2+1) %>%
					arrange(C2) %>%
					mutate(C = 1:n()) %>%
					select(name, C) %>%
					rename(`Cell type category` = name)					
				# tree %>%
				# 	ToDataFrameTree("name", "level", sprintf("C%s", 2)) %>%
				# 	arrange(C2) %>%
				# 	drop_na() %>%
				# 	select(name, C2) %>%
				# 	rename(`Cell type category` = name) %>%
				# 	rename(C = C2)
			)
		
		alpha_2 = alpha_2 %>%
			separate(.variable, c("par", "node"), remove = F) %>%
			left_join(
				fit2 %>% tidybayes::gather_draws(`prop_[a-z]_rng`[Q, C], regex = T) %>% ungroup() %>% mutate(.variable = gsub("_rng", "", .variable)) %>% separate(.variable, c("par", "node"), remove = F)  %>% select(-par) %>% nest(rng = -c(node, C)) 
			)
		
		
		internals$alpha = internals$alpha  %>% bind_rows(alpha_2 %>% mutate(level = 2))
		
		
	}
	
	internals$prop = bind_rows(internals$prop , prop_2) 
	internals$fit = internals$fit %>% c(list(fit2))
	internals$df = internals$df %>% c(list(df2))
	internals$prop_posterior[[2]] =  fit_prop_parsed %>% prop_to_list %>% `[[` (sprintf("prop_%s", "a"))  
	internals$draws = internals$draws %>% c(list(draws_2))

	internals
	
}



gather_graws_multi = function(.fit, pars){
	
	map_dfr(
		pars,
		~.fit %>% tidybayes::gather_draws(.x[Q, C], regex = T)
	)
	tidybayes::gather_draws(`prop_[a-z]`[Q, C], regex = T) %>%
		bind_rows(	tidybayes::gather_draws(`prop_[a-z]`[Q, C], regex = T))
}
#------------------------------------#

run_lv_3 = function(internals,
										shards,
										level,
										full_bayesian,
										approximate_posterior,
										iterations = iterations,
										sampling_iterations = sampling_iterations	,
										do_regression = do_regression,
										family = family,
										.formula = .formula, model = stanmodels$ARMET_tc_fix_hierarchical){
	# browser()
	res3 = run_model(
		internals$reference_filtered,
		internals$mix,
		shards,
		level,
		full_bayesian,
		approximate_posterior,
		internals$prop_posterior,
		draws_to_exposure(internals$fit[[1]]),
		iterations = iterations,
		sampling_iterations = sampling_iterations	,
		X = internals$X,
		do_regression = do_regression,
		family = family,
		cens = internals$cens,
		tree_properties = internals$tree_properties,
		Q = internals$Q,
		model = model
	)
	
	df3 = res3[[1]]
	fit3 = res3[[2]]
	
	fit_prop_parsed = fit3 %>%
		######## ALTERED WITH TREE
		tidybayes::gather_draws(prop_b[Q, C], prop_c[Q, C], prop_d[Q, C], prop_e[Q, C], prop_f[Q, C])
	
	draws_3 =
		internals$draws[[2]] %>%
		left_join(
			fit_prop_parsed %>%
			drop_na %>%
				ungroup() %>%
				left_join(
					######## ALTERED WITH TREE
					tibble(
						.variable = c("prop_b", "prop_c", "prop_d", "prop_e", "prop_f"),
						C2 = internals$tree_properties$parents_lv3
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
	
	if (do_regression && length(parse_formula(.formula)$covariates) >0 ){
		alpha_3 =
			fit3 %>%
			tidybayes::gather_draws(alpha_b[A, C], alpha_c[A, C], alpha_d[A, C], alpha_e[A, C], alpha_f[A, C]) %>%
			ungroup() %>%
			
			drop_na() %>%
			arrange(.variable) %>%
			
			# rebuild the last component sum-to-zero
			ifelse_pipe(family == "dirichlet" | 1, ~ .x %>% rebuild_last_component_sum_to_zero) %>%
			
			# Calculate relative 0 because of dirichlet relativity
			ifelse_pipe(family == "dirichlet" | 1, ~ .x %>% get_relative_zero, ~ .x %>% mutate(zero = 0)) %>%
			
			arrange(.chain, .iteration, .draw,     A) %>%
			
			nest(draws = -c(C, .variable, zero)) %>%
			
			# Attach convergence information
			left_join(
				rstan::summary(fit3, pars=sprintf("alpha_%s", c("b", "c", "d", "e", "f")))$summary %>% 
					as_tibble(rownames=".variable") %>% 
					separate(.variable, c(".variable", "A", "C"), sep="\\[|\\]|,", extra = "drop") %>% 
					mutate(A = as.integer(A), C = as.integer(C)) %>% filter(A == 2) %>% 
					select(.variable, C, Rhat),
				by = c(".variable", "C")
			) %>%
			
			# FOR HIERARCHICAL
			mutate(C = 1:n()) %>%
			
			left_join(
				tree %>%
					ToDataFrameTree("name", "level", sprintf("C%s", 3)) %>%
					filter(level == 3+1) %>%
					arrange(C3) %>%
					mutate(C = 1:n()) %>%
					select(name, C) %>%
					rename(`Cell type category` = name)
				
				# tree %>%
				# 	ToDataFrameTree("name", "level", sprintf("C%s", 3)) %>%
				# 	arrange(C3) %>%
				# 	drop_na() %>%
				# 	select(name, C3) %>%
				# 	rename(`Cell type category` = name) %>%
				# 	rename(C = C3)
			)
		
		alpha_3 = alpha_3 %>%
			separate(.variable, c("par", "node"), remove = F) %>%
			left_join(
				fit3 %>% 
					tidybayes::gather_draws(prop_b_rng[Q, C], prop_c_rng[Q, C], prop_d_rng[Q, C], prop_e_rng[Q, C], prop_f_rng[Q, C]) %>% 
					ungroup() %>%
					mutate(.variable = gsub("_rng", "", .variable)) %>% 
					separate(.variable, c("par", "node"), remove = F)  %>% 
					select(-par) %>% 
					drop_na() %>%
					nest(rng = -c(node, C)) %>%
					mutate(C = 1:n()) 
			)
		
		
		internals$alpha = 	internals$alpha  %>% bind_rows(alpha_3 %>% mutate(level = 3))
		
	}
	
	internals$prop = bind_rows(internals$prop , prop_3) 
	internals$fit = internals$fit %>% c(list(fit3))
	internals$df = internals$df %>% c(list(df3))
	internals$prop_posterior[3: (2+length(c("b", "c", "d", "e", "f")))] = fit_prop_parsed %>% prop_to_list 
	internals$draws = internals$draws %>% c(list(draws_3))

	internals
	

}

#------------------------------------#

run_lv_4 = function(internals,
										shards,
										level,
										full_bayesian,
										approximate_posterior,
										iterations = iterations,
										sampling_iterations = sampling_iterations	,
										do_regression = do_regression,
										family = family,
										.formula = .formula, model = stanmodels$ARMET_tc_fix_hierarchical){
	# browser()
	res = run_model(
		internals$reference_filtered,
		internals$mix,
		shards,
		level,
		full_bayesian,
		approximate_posterior,
		internals$prop_posterior,
		draws_to_exposure(internals$fit[[1]]),
		iterations = iterations,
		sampling_iterations = sampling_iterations	,
		X = internals$X,
		do_regression = do_regression,
		family = family,
		cens = internals$cens,
		tree_properties = internals$tree_properties,
		Q = internals$Q,
		model = model
	)
	
	df = res[[1]]
	fit = res[[2]]
	
	fit_prop_parsed = fit %>%
		######## ALTERED WITH TREE
		tidybayes::gather_draws(prop_g[Q, C], prop_h[Q, C], prop_i[Q, C], prop_l[Q, C], prop_m[Q, C])
		#########################
		
	draws =
		internals$draws[[level-1]] %>%
		left_join(
			fit_prop_parsed %>%
			drop_na %>%
				ungroup() %>%
				left_join(
					######## ALTERED WITH TREE
					tibble(
						.variable = c("prop_g", "prop_h", "prop_i", "prop_l", "prop_m"),
						C3 = internals$tree_properties$parents_lv4
					)
					#########################
				) %>%
				select(-.variable) %>%
				rename(.value4 = .value, C4 = C),
			by = c(".chain", ".iteration", ".draw", "Q",  "C3")
		) %>%
		group_by(.chain, .iteration, .draw, Q) %>%
		arrange(C1, C2, C3, C4) %>%
		mutate(
			C4 = tree$Get(sprintf("C%s", level)) %>% na.omit,
			`Cell type category` = tree$Get(sprintf("C%s", level)) %>% na.omit %>% names
		) %>%
		ungroup() %>%
		mutate(.value_relative = .value4) %>%
		mutate(.value4 = ifelse(.value4 %>% is.na, .value2, .value3 * .value4))
	
	prop =
		draws %>%
		select(.chain,
					 .iteration,
					 .draw,
					 Q,
					 C4 ,
					 `Cell type category`,
					 .value4,
					 .value_relative) %>%
		rename(C = C4, .value = .value4) %>%
		mutate(.variable = sprintf("prop_%s", level)) %>%
		mutate(level := !!level) %>%
		
		# add sample annotation
		left_join(df %>% distinct(Q, sample), by = "Q")	%>%
		
		# If MCMC is used check divergencies as well
		ifelse_pipe(
			!approximate_posterior,
			~ .x %>% parse_summary_check_divergence(),
			~ .x %>% parse_summary() %>% rename(.value = mean)
		) %>%
		
		left_join(df %>% distinct(Q, sample))
	
	if (do_regression && length(parse_formula(.formula)$covariates) >0 ){
		alpha_4 =
			fit %>%
			tidybayes::gather_draws(alpha_g[A, C], alpha_h[A, C], alpha_i[A, C], alpha_l[A, C], alpha_m[A, C]) %>%
			ungroup() %>%
			
			# rebuild the last component sum-to-zero
			ifelse_pipe(family == "dirichlet" | 1, ~ .x %>% rebuild_last_component_sum_to_zero) %>%
			
			# Calculate relative 0 because of dirichlet relativity
			ifelse_pipe(family == "dirichlet" | 1, ~ .x %>% get_relative_zero, ~ .x %>% mutate(zero = 0)) %>%
			
			arrange(.chain, .iteration, .draw,     A) %>%
			
			nest(draws = -c(C, .variable, zero)) %>%
			
			# Attach convergence information
			left_join(
				rstan::summary(fit, pars=sprintf("alpha_%s", c("g", "h", "i", "l", "m")))$summary %>% 
					as_tibble(rownames=".variable") %>% 
					separate(.variable, c(".variable", "A", "C"), sep="\\[|\\]|,", extra = "drop") %>% 
					mutate(A = as.integer(A), C = as.integer(C)) %>% filter(A == 2) %>% 
					select(.variable, C, Rhat),
				by = c(".variable", "C")
			) %>%
			
			# FOR HIERARCHICAL
			mutate(C = 1:n()) %>%
			
			left_join(
				
				tree %>%
					ToDataFrameTree("name", "level", sprintf("C%s", 4)) %>%
					filter(level == 4+1) %>%
					arrange(C4) %>%
					mutate(C = 1:n()) %>%
					select(name, C) %>%
					rename(`Cell type category` = name)
				
				# tree %>%
				# 	ToDataFrameTree("name", "level", sprintf("C%s", 4)) %>%
				# 	arrange(C4) %>%
				# 	drop_na() %>%
				# 	select(name, C4) %>%
				# 	rename(`Cell type category` = name) %>%
				# 	rename(C = C4)
			)
		
		alpha_4 = alpha_4 %>%
			separate(.variable, c("par", "node"), remove = F) %>%
			left_join(
				fit %>% 
					tidybayes::gather_draws(prop_g_rng[Q, C], prop_h_rng[Q, C], prop_i_rng[Q, C], prop_l_rng[Q, C], prop_m_rng[Q, C]) %>% 
					ungroup() %>% mutate(.variable = gsub("_rng", "", .variable)) %>% 
					separate(.variable, c("par", "node"), remove = F)  %>% 
					select(-par) %>% drop_na() %>% nest(rng = -c(node, C)) %>%
					mutate(C = 1:n()) 
			)
		
		
		internals$alpha = 	internals$alpha  %>% bind_rows(alpha_4 %>% mutate(level = 4))
		
		
		
	}
	
	internals$prop = bind_rows(internals$prop , prop) 
	internals$fit = internals$fit %>% c(list(fit))
	internals$df = internals$df %>% c(list(df))
	internals$draws = internals$draws %>% c(list(draws))
	
	
	internals

	
}
