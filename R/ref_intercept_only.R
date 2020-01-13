
#' ref_intercept_only
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
ref_intercept_only = function(
	reference,
	level,
	cores = 8
){

	# Former parameters
	X = matrix(rep(1, 1))
	do_regression = F
	full_bayesian = T
	approximate_posterior = F
	omit_regression =                   T
	save_fit =                          T
	seed =                              NULL
	iterations = 250
	sampling_iterations = 100
	levels = 1:4


	mix = reference %>% group_by(symbol) %>%
		summarise(count = `read count` %>% median(na.rm = T) %>% as.integer) %>%
		mutate(sample = "dummy") %>%
		spread(symbol, count)

	# Add fake lambda, sigma
	reference =
		reference %>%
		mutate(lambda = 1, sigma_raw = 1)

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

	tree = create_tree_object(reference)

	tree = 	data.tree::Clone(tree) %>%	{
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
		reference %>%
		select(level, sample, symbol, `Cell type category`, `read count`, lambda, sigma_raw, `house keeping`) %>%

		# Bug after I deleted FANTOM5 I have to rerun infer NB. Some genes are not in all cell types anynore
		# Other bug
		filter(`Cell type category` %>% is.na %>% `!`) %>%

		# Check if this is still important
		group_by(level) %>%
		do(
			(.) %>% inner_join(
				(.) %>%
					distinct(symbol, `Cell type category`) %>%
					count(symbol) %>%
					filter(n == max(n))
			)
		) %>%
		ungroup() %>%


		# left_join(n_markers, by=c("ct1", "ct2")) %>%
		# filter_reference(mix) %>%
		# select(-ct1, -ct2, -rank, -`n markers`) %>%
		# distinct %>%

		# Select cell types in hierarchy
		inner_join(
			tree %>%
				data.tree::ToDataFrameTree("Cell type category", "C", "C1", "C2", "C3") %>%
				as_tibble %>%
				select(-1)

		)

	# # CAP THE SIGMA TO AVOID OVERFITTING
	# reference_filtered = reference_filtered %>%
	# 	#mutate(sigma_raw_capped = ifelse(sigma_raw > sigma_raw_minimum, sigma_raw, sigma_raw_minimum)) %>%
	# 	mutate(sigma_raw_capped = ifelse(sigma_raw > sigma_raw_regressed, sigma_raw, sigma_raw_regressed)) %>%
	# 	mutate(sigma_raw = sigma_raw_capped)


	prop_posterior = get_null_prop_posterior(ct_in_nodes)

	######################################

	res1 = run_model(
		reference_filtered,
		mix,
		shards,
		level,
		T,
		approximate_posterior,
		prop_posterior,
		iterations = iterations,
		sampling_iterations = sampling_iterations,
		X = X, do_regression = do_regression
	)

	res1[[1]] %>% filter(!query) %>% distinct(symbol, `Cell type category`, G) %>%
		left_join(
			res1[[2]] %>% rstan::summary(c("lambda_log", "sigma_inv_log")) %$% summary %>%
		as_tibble(rownames="par") %>%
		separate(par, c("par", "G"), sep="\\[|\\]", extra = "drop") %>%
		mutate(G = G %>% as.integer) %>%
		select(par, G, mean) %>%
		spread(par, mean),
		by = c("G")
	)
}







