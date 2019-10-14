
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
	approximate_posterior = F,
	verbose =                           F,
	omit_regression =                   F,
	save_fit =                          F,
	seed =                              NULL,
	cores = 14,
	iterations = 300,
	sampling_iterations = 100,
	levels = 1:4,
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
			#aspect.ratio=1,
			axis.text.x = element_text(angle = 30, hjust = 1),
			strip.background = element_blank(),
			axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
			axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
		)

	shards = cores #* 2
	shards = cores = 8
	shards = shards
	is_level_in = shards %>% `>` (0) %>% as.integer

	# Global properties - derived by previous analyses of the whole reference dataset
	sigma_intercept = 1.3420415
	sigma_slope = -0.3386389
	sigma_sigma = 1.1720851
	lambda_mu_mu = 5.612671
	lambda_sigma = 7.131593

	# Non centered
	lambda_mu_prior = c(6.2, 1)
	lambda_sigma_prior =  c( log(3.3) , 1)
	lambda_skew_prior =  c( -2.7, 1)
	sigma_intercept_prior = c( 1.9 , 0.1)

	# Set up tree structure
	levels_in_the_tree = 1:4

	tree = 	data.tree::Clone(tree) %>%	{
		# Filter selected levels
		data.tree::Prune(., function(x) x$level <= max(levels_in_the_tree) + 1)
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

	Q = mix %>% nrow

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

		)
	#%>%
	# # Decrease the number of house keeping used
	# anti_join({
	# 	mdf = (.) %>%
	# 		distinct(symbol, `house keeping`) %>%
	# 		filter(`house keeping`)
	#
	# 	withr::with_seed(	123, 	sample_frac(mdf, 0.8)) %>%
	# 		distinct(symbol)
	# }) %>%

	prop_posterior = get_null_prop_posterior(ct_in_nodes)


	######################################

	res1 = run_model(	reference_filtered, mix, shards,	1,	full_bayesian, approximate_posterior, prop_posterior	)

	df1 = res1[[1]]
	fit1 = res1[[2]]

	prop_posterior[[1]] = fit1 %>% draws_to_alphas("prop_1") %>% `[[` (1)

	prop1 = get_prop(fit1, approximate_posterior, df1)

	draws_1 =
		fit1 %>%
		tidybayes::gather_draws(prop_1[Q, C]) %>%
		ungroup() %>%
		select(-.variable)

	######################################

	res2 = run_model(	reference_filtered, mix, shards,	2,	full_bayesian, approximate_posterior, prop_posterior, draws_to_exposure(fit1)	)

	df2 = res2[[1]]
	fit2 = res2[[2]]

	fit2_to_prop2 = function(){

	}

	prop_2 =
		draws_1 %>%
		left_join(
			fit2 %>%
			tidybayes::gather_draws(prop_a[Q, C]) %>%
			ungroup() %>%
			select(-.variable) %>%
			rename(C2 = C, .value2 = .value) %>%
			mutate(C = parents_lv2[1]),
		by = c(".chain", ".iteration", ".draw", "Q",  "C")
		) %>%
		mutate(
			.value2 = ifelse(.value2 %>% is.na, .value, .value * .value2),
			C2 = ifelse(C2 %>% is.na, C, C + C2 -1)
		) %>%
		select(-C, -.value) %>%
		rename(C = C2, .value = .value2) %>%
		mutate(.variable = "prop_2") %>%
		group_by(.variable,  Q,  C) %>%
		tidybayes::median_qi()

	prop_posterior[[2]] = fit2 %>% draws_to_alphas(sprintf("prop_%s", "a")) %>% `[[` (1)


	###########################################

	browser()

	res3 = run_model(	reference_filtered, mix, shards,	3,	full_bayesian, approximate_posterior, prop_posterior, draws_to_exposure(fit2)	)

	df3 = res3[[1]]
	fit3 = res3[[2]]


	# Produce results

	# Return
	list(

		# Matrix of proportions
		proportions =	prop,

		# Return the input itself
		input = input,

		# Return the fitted object
		fit = fit,

		# # Return data source
		# data_source = y_source,
		signatures = reference_filtered,
		signatures1 = df1
	)

}



