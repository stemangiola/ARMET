#' This is a generalisation of ifelse that acceots an object and return an objects
#'
#' @import dplyr
#' @import tidyr
#'
#' @param input.df A tibble
#' @param condition A boolean
#' @return A tibble
ifelse_pipe = function(.x, .p, .f1, .f2 = NULL) {
	switch(.p %>% `!` %>% sum(1),
				 as_mapper(.f1)(.x),
				 if (.f2 %>% is.null %>% `!`)
				 	as_mapper(.f2)(.x)
				 else
				 	.x)

}


#' format_for_MPI
#'
#' @description Format reference data frame for MPI
format_for_MPI = function(df, shards){
	df %>%

		left_join(
			(.) %>%
				distinct(G) %>%
				arrange(G) %>%
				mutate( idx_MPI = head( rep(1:shards, (.) %>% nrow %>% `/` (shards) %>% ceiling ), n=(.) %>% nrow) )
		) %>%
		arrange(idx_MPI, G) %>%

		# Decide start - end location
		group_by(idx_MPI) %>%
		do(
			(.) %>%
				left_join(
					(.) %>%
						distinct(sample, G) %>%
						arrange(G) %>%
						count(G) %>%
						mutate(end = cumsum(n)) %>%
						mutate(start = c(1, .$end %>% rev() %>% `[` (-1) %>% rev %>% `+` (1)))
				)
		) %>%
		ungroup() %>%

		# Add ct_symbol MPI rows indexes - otherwise spread below gives error
		left_join(
			(.) %>%
				group_by(idx_MPI) %>%
				distinct(G) %>%
				arrange(G) %>%
				mutate(`symbol MPI row` = 1:n()) %>%
				ungroup
		) %>%

		# Add counts MPI rows indexes
		group_by(idx_MPI) %>%
		arrange(G) %>%
		mutate(`read count MPI row` = 1:n()) %>%
		ungroup

}

#' add_partition
#'
#' @description Add partition column dto data frame
add_partition = function(df.input, partition_by, n_partitions){
	df.input %>%
		left_join(
			(.) %>%
				select(!!partition_by) %>%
				distinct %>%
				mutate(
					partition = 1:n() %>%
						divide_by(length((.))) %>%
						#	multiply_by(min(n_partitions, df.input %>% distinct(symbol) %>% nrow)) %>%
						multiply_by(n_partitions) %>%
						ceiling
				)
		)
}

#' get_MPI_deconv
#'
#' @description Get data format for MPI deconvolution part
get_MPI_deconv = function(y_source, shards, my_level, tree){
	y_MPI_source =
		y_source %>%
		distinct(level, Q, S, symbol, G, GM, `Cell type category`, `read count`) %>%
		filter(level==my_level)  %>%

		# Add universal cell type rank
		left_join(
			tree %>%
				ToDataFrameTree("name", sprintf("C%s", my_level)) %>%
				select(-1) %>%
				setNames(c("Cell type category" , "ct_rank"))
		) %>%

		# Arrange very important for consistency
		arrange(Q, symbol, ct_rank) %>%
		add_partition("symbol", shards) %>%
		group_by(partition) %>%
		left_join( (.) %>% distinct(symbol) %>% mutate(MPI_row = 1:n())) %>%
		ungroup()

	list(
		y_MPI_source = y_MPI_source,

		y_MPI_symbol_per_shard =
			y_MPI_source %>%
			distinct(symbol, partition) %>%
			count(partition) %>%
			spread(partition, n) %>%
			as_vector %>% array,

		y_MPI_idx_symbol =
			y_MPI_source %>%
			distinct(MPI_row, GM, partition, Q) %>%
			spread(partition, GM) %>%
			select(-Q) %>% distinct %>%
			replace(is.na(.), 0 %>% as.integer) %>%
			as_matrix(rownames = "MPI_row") %>%
			t,

		y_MPI_G_per_shard =
			y_MPI_source %>%
			distinct(symbol, `Cell type category`, partition) %>%
			count(partition) %>%
			spread(partition, n) %>%
			as_vector %>% array,

		y_MPI_idx =
			y_MPI_source %>%
			distinct(partition, symbol, G, `Cell type category`, MPI_row) %>%
			select(-symbol) %>%
			spread(partition, G) %>%
			arrange(MPI_row, `Cell type category`) %>%
			select(-MPI_row,-`Cell type category`) %>%
			replace(is.na(.), 0 %>% as.integer) %>%
			as_matrix %>%
			t,

		y_idx =
			y_MPI_source %>% distinct(symbol, G, `Cell type category`) %>% arrange(symbol, `Cell type category`) %>% pull(G),

		y_MPI_N_per_shard =
			y_MPI_source %>%
			distinct(MPI_row, `read count`, partition, Q) %>%
			count(partition) %>%
			spread(partition, n) %>%
			as_vector %>% array,

		y_MPI_count =
			y_MPI_source %>%
			distinct(MPI_row, `read count`, partition, Q) %>%
			spread(partition, `read count`) %>%
			arrange(MPI_row, Q) %>%
			select(-MPI_row,-Q) %>%
			replace(is.na(.), 0 %>% as.integer) %>%
			as_matrix %>%
			t
	)
}

plot_differences_in_lambda = function(){

	# Plot differences in lambda
	(fit %>%
	 	tidybayes::gather_draws(lambda_log[G]) %>%
	 	tidybayes::median_qi() %>%
	 	left_join(
	 		counts_baseline %>%
	 			distinct(`symbol`, G, `Cell type category`)
	 	) %>%
	 	left_join(
	 		reference_filtered %>%
	 			distinct(symbol, lambda, `Cell type category`) %>%
	 			rename(`lambda_log` = lambda)
	 	) %>%
	 	ggplot(aes(x=lambda_log, y=.value, label=G)) + geom_point() + geom_abline(intercept = 0, slope = 1, color="red") + my_theme
	)  %>% plotly::ggplotly()

	#
	(
		fit %>%
			extract(pars=c("lambda_mu", "lambda_sigma", "exposure_rate",  "lambda_log", "sigma_raw", "prop")) %>%
			as.data.frame %>% as_tibble() %>%
			mutate(chain = rep(1:3, 100) %>% sort %>% as.factor ) %>%
			select(chain, everything()) %>% gather(par, draw, -chain) %>%
			group_by(chain, par) %>%
			summarise(d = draw %>% median) %>%
			ggplot(aes(y=d, x=par, color=chain)) + geom_point()
	) %>% plotly::ggplotly()

}

get_overlap_descriptive_stats = function(mix_tbl, ref_tbl){

	writeLines(
		sprintf(
			"%s house keeping genes are missing from the input mixture",
			ref_tbl %>% filter(`house keeping`) %>% distinct(symbol) %>% anti_join( mix_tbl %>% distinct(symbol) ) %>% nrow
		)
	)

	writeLines(
		sprintf(
			"%s marker genes are missing from the input mixture",
			ref_tbl %>% filter(!`house keeping`) %>% distinct(symbol) %>% anti_join( mix_tbl %>% distinct(symbol) ) %>% nrow
		)
	)

	ref_tbl %>%
		filter(!`house keeping`) %>% distinct(symbol, ct1, ct2) %>% anti_join( mix_tbl %>% distinct(symbol)) %>%
		count(ct1, ct2) %>%
		rename(`missing markers` = n) %>%
		print(n=999)

}

plot_counts_inferred_sum = function(fit_obj, samples = NULL){

	fit_obj %$% fit %>%
		summary(par=c("nb_sum")) %$% summary %>%
		as_tibble(rownames="par") %>% select(par, `2.5%`, `50%`, `97.5%`) %>%
		separate(par, c(".variable", "Q", "GM"), sep="\\[|,|\\]") %>%
		mutate(Q = Q %>% as.integer, GM = GM %>% as.integer) %>%
		left_join(fit_obj %$% data_source %>% distinct(Q, GM, symbol, `read count`, sample)) %>%

		# Select samples
		{
			if(samples %>% is.null %>% `!`) (.) %>% filter(sample %in% !samples)
			else (.)
		} %>%

		# Check if inside
		rowwise %>%
		mutate(inside = between(`read count`, `2.5%`, `97.5%`)) %>%
		ungroup %>%
		ggplot(aes(x=`read count` + 1, y=`50%` + 1, color=inside, label = symbol)) +
		geom_point(alpha=0.5)  +
		geom_abline(slope = 1, intercept = 0) +
		geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), alpha=0.5) +
		facet_wrap(~sample) +
		scale_y_log10() +
		scale_x_log10()

}

choose_chains_majority_roule = function(fit_parsed){

	fit_parsed %>%
		inner_join(

			# Calculate modes
			fit_parsed %>%
				select(.chain, .value) %>%
				{
					main_cluster =
						(.) %>%
						pull(.value) %>%
						kmeans(centers = 2) %>% {
							bind_cols(
								(.) %$% centers %>% as_tibble(.name_repair = "minimal") %>% setNames("center") ,
								(.) %$% size %>% as_tibble(.name_repair = "minimal")  %>% setNames("size")
							)
						} %>%
						arrange(size %>% desc) %>%
						slice(1)

					(.) %>%
						group_by(.chain) %>%
						summarise(.lower_chain = quantile(.value, probs = c(0.025)), .upper_chain = quantile(.value, probs = c(0.975))) %>%
						ungroup %>%
						mutate(center = main_cluster %>% pull(center))
				} %>%

				# Filter cains
				rowwise() %>%
				filter(between(center, .lower_chain, .upper_chain)) %>%
				ungroup %>%
				distinct(.chain)

		)
}

filter_reference = function(reference, mix){
	reference %>%
		inner_join(mix %>% slice(1) %>%	gather(`symbol`, `read count`, -sample) %>% distinct(symbol)) %>%
		{
			bind_rows(

				# Get markers based on common genes with mix
				(.) %>%
					filter(!`house keeping`) %>%
					group_by(ct1, ct2) %>%
					do({

						n_markers = (.) %>% slice(1) %>% pull(`n markers`)
						(.) %>%
							inner_join(
								(.) %>%
									distinct(symbol, rank) %>%
									arrange(rank) %>%
									slice(1:n_markers)
							)

					}) %>%
					ungroup() %>%

					# Filter markers in upper levels that have been selected for lower levels
					inner_join(
						(.) %>%
							distinct(symbol, level) %>%
							group_by(symbol) %>%
							arrange(level %>% desc) %>%
							slice(1) %>%
							ungroup
					) %>%

					# Print number of markewrs per comparison
					{
						(.) %>% distinct(symbol, ct1, ct2) %>% count(ct1, ct2) %>% print(n=99)
						(.)
					},

				# Get house keeping genes
				(.) %>% filter(`house keeping`)
			)
		} %>%
		select(-ct1, -ct2, -rank, -`n markers`) %>%	distinct
}

get_idx_level = function(tree, my_level){
	left_join(
		tree %>% ToDataFrameTree("name") %>% as_tibble,
		Clone(tree) %>%
			{ Prune(., function(x) x$level <= my_level + 1); . } %>%
			ToDataFrameTree("level", "C", "isLeaf", "name") %>%
			as_tibble %>%
			filter(isLeaf) %>%
			left_join(

				tree %>%
					ToDataFrameTree("name", "level", "isLeaf") %>%
					as_tibble %>%
					filter(level <= my_level & isLeaf) %>%
					mutate(isAncestorLeaf = !isLeaf) %>%
					select(name, isAncestorLeaf)

			) %>%
			arrange(isAncestorLeaf) %>%
			mutate(my_C = 1:n()) %>%
			select(name, my_C)
	) %>%
		pull(my_C)
}

parse_summary = function(fit){
	fit %>%
		rstan::summary() %$% summary %>%
		as_tibble(rownames=".variable") %>%
		filter(grepl("prop", .variable)) %>%
		separate(.variable, c(".variable", "Q", "C"), sep="[\\[,\\]]", extra="drop") %>%
		mutate(C = C %>% as.integer, Q = Q %>% as.integer)
}

parse_summary_check_divergence = function(fit){
	fit %>%
		tidybayes::gather_draws(prop_1[Q, C], prop_2[Q, C], prop_3[Q, C]) %>%
		filter(.variable %in% c("prop_1", "prop_2", "prop_3")) %>%

		# If not converged choose the majority chains
		mutate(	converged = diptest::dip.test(`.value`) %$%	`p.value` > 0.05) %>%

		# If some proportions have not converged chose the most populated one
		do(
			(.) %>%
				ifelse_pipe(
					(.) %>% distinct(converged) %>% pull(1) %>% `!`,
					~ .x %>% choose_chains_majority_roule
				)
		) %>%

		# Anonymous function - add summary fit to converged label
		# input: tibble
		# output: tibble
		{
			left_join(
				(.) %>% select(-converged) %>% tidybayes::median_qi(),
				(.) %>% distinct(converged)
			)
		} %>%
		ungroup()
}


library(data.tree)
tree =
	yaml:: yaml.load_file("~/PhD/deconvolution/ARMET/data/tree.yaml") %>%
	data.tree::as.Node() %>%

	{

		# Sort tree by name
		Sort(., "name")

		# Add C indexes
		.$Set(
			C =
				tibble( name = .$Get('name'), level = .$Get('level')) %>%
				left_join( (.) %>% arrange(level, name) %>%	 	mutate(C = 0:(n()-1))  )	%>%
				pull(C)
		)
		.$Set(	C1 = get_idx_level(.,1)	)
		.$Set(	C2 = get_idx_level(.,2)	)
		.$Set(	C3 = get_idx_level(.,3)	)
		#		if(max(levels)>1) for(l in 2:max(levels)) { my_c = sprintf("C%s", l); .$Set(	my_c = get_idx_level(.,2)	); . }

		# Set Cell type category label
		.$Set("Cell type category" = .$Get("name"))

	}
library(magrittr)
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
#' @param mix A matrix
#' @param my_design A matrix
#' @param cov_to_test A character string
#' @param fully_bayesian A boolean
#' @param is_mix_microarray A boolean
#' @param ct_to_omit A character string
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
	reference,
	full_bayesian = 0,
	ct_to_omit =                        c("t_CD4_naive", "adipocyte"),
	verbose =                           F,
	omit_regression =                   F,
	save_fit =                          F,
	seed =                              NULL,
	cores = 14,
	iterations = 300,
	sampling_iterations = 100,
	levels = 1:3,
	full_bayes = T
){

	# full_bayesian = 0
	# ct_to_omit =                        c("t_CD4_naive", "adipocyte")
	# verbose =                           F
	# omit_regression =                   F
	# save_fit =                          F
	# seed =                              NULL
	# cores = 14
	#levels = 1:2

	source("https://gist.githubusercontent.com/stemangiola/dd3573be22492fc03856cd2c53a755a9/raw/e4ec6a2348efc2f62b88f10b12e70f4c6273a10a/tidy_extensions.R")
	source("https://gist.githubusercontent.com/stemangiola/90a528038b8c52b21f9cfa6bb186d583/raw/4a5798857362d946bd3029188b1cc9eb9b625456/transcription_tool_kit.R")
	source("https://gist.githubusercontent.com/stemangiola/9d2ba5d599b7ac80404c753cdee04a01/raw/ad571ea2bbc3f13441a7d845b5ae8ed67a45d8ec/tidy_data_tree.R")

	library(tidyverse)
	library(foreach)
	library(rstan)

	input = c(as.list(environment()))
	shards = cores #* 2

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

	#########################################
	# Set up tree structure
	#########################################

	tree = 	tree %>%	{
		# Filter selected levels
			Prune(., function(x) x$level <= max(levels) + 1)
			.
		}

	ct_in_nodes = tree %>% ToDataFrameTree("name", "level", "C", "count", "isLeaf") %>% as_tibble %>% arrange(level, C) %>% filter(!isLeaf) %>% pull(count)
	ct_in_levels = foreach(l=levels+1, .combine = c) %do% {

	Clone(tree) %>%
		{ Prune(., function(x) x$level <= l);	. 	} %>%
		ToDataFrameTree("name", "level", "C", "count", "isLeaf") %>%
		as_tibble %>%
		arrange(level, C) %>%
		filter(isLeaf) %>%
		nrow
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

	# Print overlap descriptive stats
	#get_overlap_descriptive_stats(mix %>% slice(1) %>% gather(`symbol`, `read count`, -sample), reference)

	#########################################
	# Prepare data frames -
	# For Q query first
	# For G house keeing first
	# For GM level 1 first
	#########################################

	reference_filtered =
		filter_reference(reference, mix) %>%

		# Select cell types in hierarchy
		inner_join(
			tree %>%
				ToDataFrameTree("Cell type category", "C") %>%
				as_tibble %>%
				select(-1)

		)

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
				slice(2:n()) %>%
				ungroup()
		) %>%

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

	#########################################
	# For  reference MPI inference
	#########################################

	counts_baseline =
		df %>%

		# Eliminate the query part, not the house keeping of the query
		filter(!`query` | `house keeping`)  %>%

		# If full Bayesian false just keep house keeping
		{
			if(!full_bayesian) (.) %>% filter(`house keeping` & `query`)
			else (.)
		} %>%

		format_for_MPI(shards)

	S = counts_baseline %>% distinct(sample) %>% nrow()
	N = counts_baseline %>% distinct(idx_MPI, `read count`, `read count MPI row`) %>%  count(idx_MPI) %>% summarise(max(n)) %>% pull(1)
	M = counts_baseline %>% distinct(start, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% max
	G_per_shard = counts_baseline %>% distinct(ct_symbol, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% as.array
	n_shards = min(shards, counts_baseline %>% distinct(idx_MPI) %>% nrow)
	G_per_shard_idx = c(0, counts_baseline %>% distinct(ct_symbol, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% cumsum)

	counts =
		counts_baseline %>%
		distinct(idx_MPI, `read count`, `read count MPI row`)  %>%
		spread(idx_MPI,  `read count`) %>%
		select(-`read count MPI row`) %>%
		replace(is.na(.), 0 %>% as.integer) %>%
		as_matrix() %>% t

	sample_idx =
		counts_baseline %>%
		distinct(idx_MPI, S, `read count MPI row`)  %>%
		spread(idx_MPI, S) %>%
		select(-`read count MPI row`) %>%
		replace(is.na(.), 0 %>% as.integer) %>%
		as_matrix() %>% t

	symbol_end =
		counts_baseline %>%
		distinct(idx_MPI, end, `symbol MPI row`)  %>%
		spread(idx_MPI, end) %>%
		bind_rows( (.) %>% head(n=1) %>%  mutate_all(function(x) {0}) ) %>%
		arrange(`symbol MPI row`) %>%
		select(-`symbol MPI row`) %>%
		replace(is.na(.), 0 %>% as.integer) %>%
		as_matrix() %>% t

	G_ind =
		counts_baseline %>%
		distinct(idx_MPI, G, `symbol MPI row`)  %>%
		spread(idx_MPI, G) %>%
		arrange(`symbol MPI row`) %>%
		select(-`symbol MPI row`) %>%
		replace(is.na(.), 0 %>% as.integer) %>%
		as_matrix() %>% t

	#######################################
	# For deconvolution
	#######################################

	n_house_keeping = df %>% filter(!`query` & `house keeping`) %>% distinct(G, `house keeping`) %>% nrow

	y_source =
		df %>%
		filter(`query` & !`house keeping`) %>%
		select(S, Q, `symbol`, `read count`, GM, sample) %>%
		left_join(	df %>% filter(!query) %>% distinct(`symbol`, G, `Cell type category`, level, lambda, sigma_raw, GM, C) ) %>%
		arrange(C, Q, symbol) %>%
		mutate(`Cell type category` = factor(`Cell type category`, unique(`Cell type category`)))

	# Data MPI for deconvolution level 1
	y_MPI_lv1 = y_source %>% get_MPI_deconv(shards, 1, tree)
	idx_1_source = y_source %>% filter(level == 1) %>% distinct(symbol, G, `Cell type category`) %>% arrange(`Cell type category`, symbol)
	idx_1 = idx_1_source %>% pull(G)
	I1 = idx_1 %>% length
	I1_dim = c(idx_1_source %>% distinct(symbol) %>% nrow, idx_1_source %>% distinct(`Cell type category`) %>% nrow)
	y_MPI_source_lv1 = y_MPI_lv1 %$% y_MPI_source
	y_MPI_symbol_per_shard_lv1 = y_MPI_lv1 %$% y_MPI_symbol_per_shard
	y_MPI_idx_symbol_lv1 = y_MPI_lv1 %$% y_MPI_idx_symbol
	y_MPI_G_per_shard_lv1 = y_MPI_lv1 %$% y_MPI_G_per_shard
	GM_lv1 = df %>% filter(!`house keeping`) %>% filter(level==1) %>% distinct(symbol) %>% nrow()
	y_idx_lv1 =  y_MPI_lv1 %$% y_idx
	y_MPI_idx_lv1 = y_MPI_lv1 %$% y_MPI_idx
	y_MPI_N_per_shard_lv1 = y_MPI_lv1 %$% y_MPI_N_per_shard
	y_MPI_count_lv1 = y_MPI_lv1 %$% y_MPI_count

	# Data MPI for deconvolution level 2
	y_MPI_lv2 = y_source %>% get_MPI_deconv(shards, 2, tree)
	idx_2_source = y_source %>% filter(level == 2) %>% distinct(symbol, G, `Cell type category`) %>% arrange(`Cell type category`, symbol)
	idx_2 = idx_2_source %>% pull(G)
	I2 = idx_2 %>% length
	I2_dim = c(idx_2_source %>% distinct(symbol) %>% nrow, idx_2_source %>% distinct(`Cell type category`) %>% nrow)
	y_MPI_source_lv2 = y_MPI_lv2 %$% y_MPI_source
	y_MPI_symbol_per_shard_lv2 = y_MPI_lv2 %$% y_MPI_symbol_per_shard
	y_MPI_idx_symbol_lv2 = y_MPI_lv2 %$% y_MPI_idx_symbol
	y_MPI_G_per_shard_lv2 = y_MPI_lv2 %$% y_MPI_G_per_shard
	GM_lv2 = df %>% filter(!`house keeping`) %>% filter(level==2) %>% distinct(symbol) %>% nrow()
	y_idx_lv2 =  y_MPI_lv2 %$% y_idx
	y_MPI_idx_lv2 = y_MPI_lv2 %$% y_MPI_idx
	y_MPI_N_per_shard_lv2 = y_MPI_lv2 %$% y_MPI_N_per_shard
	y_MPI_count_lv2 = y_MPI_lv2 %$% y_MPI_count

	# Data MPI for deconvolution level 3
	y_MPI_lv3 = y_source %>% get_MPI_deconv(shards, 3, tree)
	idx_3_source = y_source %>% filter(level == 3) %>% distinct(symbol, G, `Cell type category`) %>% arrange(`Cell type category`, symbol)
	idx_3 = idx_3_source %>% pull(G)
	I3 = idx_3 %>% length
	I3_dim = c(idx_3_source %>% distinct(symbol) %>% nrow, idx_3_source %>% distinct(`Cell type category`) %>% nrow)
	y_MPI_source_lv3 = y_MPI_lv3 %$% y_MPI_source
	y_MPI_symbol_per_shard_lv3 = y_MPI_lv3 %$% y_MPI_symbol_per_shard
	y_MPI_idx_symbol_lv3 = y_MPI_lv3 %$% y_MPI_idx_symbol
	y_MPI_G_per_shard_lv3 = y_MPI_lv3 %$% y_MPI_G_per_shard
	GM_lv3 = df %>% filter(!`house keeping`) %>% filter(level==3) %>% distinct(symbol) %>% nrow()
	y_idx_lv3 =  y_MPI_lv3 %$% y_idx
	y_MPI_idx_lv3 = y_MPI_lv3 %$% y_MPI_idx
	y_MPI_N_per_shard_lv3 = y_MPI_lv3 %$% y_MPI_N_per_shard
	y_MPI_count_lv3 = y_MPI_lv3 %$% y_MPI_count

	Q = df %>% filter(`query`) %>% distinct(Q) %>% nrow

	#######################################
	# Merge all MPI
	#######################################

	# Reference
	counts_package =
		# Dimensions data sets
		rep(c(M, N, S), shards) %>%
		matrix(nrow = shards, byrow = T) %>%
		cbind(G_per_shard) %>%
		cbind(symbol_end) %>%
		cbind(sample_idx) %>%
		cbind(counts)

	# level 1
	lev1_package =
		# Dimensions data sets
		rep(c(ct_in_levels[1], Q, S), shards) %>%
		matrix(nrow = shards, byrow = T) %>%
		cbind(y_MPI_symbol_per_shard_lv1) %>%
		cbind(y_MPI_G_per_shard_lv1) %>%
		cbind(y_MPI_N_per_shard_lv1) %>%
		cbind(y_MPI_count_lv1)

	# level 2
	lev2_package =
		# Dimensions data sets
		rep(c(ct_in_levels[2], Q, S), shards) %>%
		matrix(nrow = shards, byrow = T) %>%
		cbind(y_MPI_symbol_per_shard_lv2) %>%
		cbind(y_MPI_G_per_shard_lv2) %>%
		cbind(y_MPI_N_per_shard_lv2) %>%
		cbind(y_MPI_count_lv2)

	# level 3
	lev3_package =
		# Dimensions data sets
		rep(c(ct_in_levels[3], Q, S), shards) %>%
		matrix(nrow = shards, byrow = T) %>%
		cbind(y_MPI_symbol_per_shard_lv3) %>%
		cbind(y_MPI_G_per_shard_lv3) %>%
		cbind(y_MPI_N_per_shard_lv3) %>%
		cbind(y_MPI_count_lv3)

	# Integrate everything
	data_package =

		# Size 3 data data sets
		rep(c(
			ncol(counts_package),
			ncol(lev1_package),
			ncol(lev2_package),
			ncol(lev3_package)
		), shards) %>%
		matrix(nrow = shards, byrow = T) %>%

		# Size parameter datasets
		cbind(
			rep(c(
				(M + 2 + S), # lambda, sigma slope intercept, exposure
				(max(y_MPI_G_per_shard_lv1) * 2 + Q + max(y_MPI_symbol_per_shard_lv1) + (Q * ct_in_levels[1])),
				(max(y_MPI_G_per_shard_lv2) * 2 + Q + max(y_MPI_symbol_per_shard_lv2) + (Q * ct_in_levels[2])),
				(max(y_MPI_G_per_shard_lv3) * 2 + Q + max(y_MPI_symbol_per_shard_lv3) + (Q * ct_in_levels[3]))
			), shards) %>%
			matrix(nrow = shards, byrow = T)
		) %>%

		# Data sets
		cbind(counts_package) %>%
		cbind(lev1_package) %>%
		cbind(lev2_package) %>%
		cbind(lev3_package)

	########################################
	########################################

	# Old data structure

	y = y_source %>% distinct(level, Q, S, symbol, `read count`) %>% arrange(level, Q, symbol) %>% select(`read count`, S) %>% as_matrix
	I = y %>% nrow


	# Pass previously infered parameters
	do_infer = full_bayesian

	lambda_log_data =
		df %>%

		# Eliminate the query part, not the house keeping of the query
		filter(!`query` | `house keeping`)  %>%

		# Ths is bcause mix lacks lambda info and produces NA in the df
		filter(!(`Cell type category` == "house_keeping" & lambda %>% is.na)) %>%

		distinct(G, lambda) %>%
		arrange(G)%>%
		pull(lambda)

	sigma_raw_data =
		df %>%

		# Eliminate the query part, not the house keeping of the query
		filter(!`query` | `house keeping`)  %>%

		# Ths is bcause mix lacks lambda info and produces NA in the df
		filter(!(`Cell type category` == "house_keeping" & sigma_raw %>% is.na)) %>%

		distinct(G, sigma_raw) %>%
		arrange(G)%>%
		pull(sigma_raw)

	########################################
	# Build better scales

	exposure_rate_shift_scale =
		df %>%
		filter(`house keeping`) %>%
		distinct(symbol, sample, `read count`, query) %>%
		drop_na %>%
		inner_join( (.) %>% distinct(sample, symbol) %>% count(symbol) %>% filter(n == max(n)), by="symbol") %>% # Eliminate genes that are missing from samples
		tidyTranscriptomics::add_normalised_counts(transcript_column = symbol) %>%
		filter(query) %>%
		mutate(l = multiplier %>% log) %>%
		summarise(shift = l %>% mean, scale = l %>% sd) %>%
		as.numeric


	intercept_shift_scale =
		counts_baseline %>%
		filter(`house keeping`) %>%
		distinct(symbol, sample, `read count`, query) %>%
		tidyTranscriptomics::add_normalised_counts(transcript_column = symbol) %>%
		mutate(
			cc = `read count normalised` %>%
				`+` (1) %>% log
		) %>%
		summarise(shift = cc %>% mean, scale = cc %>% sd) %>%
		as.numeric

	##########################################
	##########################################

	counts_linear = counts_baseline %>% arrange(G, S) %>% mutate(S = S %>% as.factor %>% as.integer) %>%  pull(`read count`)
	G_linear = counts_baseline %>% arrange(G, S) %>% mutate(S = S %>% as.factor %>% as.integer) %>% pull(G)

	S_linear = counts_baseline %>%  arrange(G, S) %>% mutate(S = S %>% as.factor %>% as.integer) %>% pull(S)
	CL = length(counts_linear)
	G = counts_baseline %>%  mutate(S = S %>% as.factor %>% as.integer)%>%  distinct(G) %>% nrow
	S = counts_baseline %>%  mutate(S = S %>% as.factor %>% as.integer)%>% distinct(S) %>% nrow

	# Deconvolution, get G only for markers of each level. Exclude house keeping
	G1_linear = counts_baseline %>% filter(!`house keeping`) %>% filter(level ==1) %>% distinct(G, GM, C) %>% arrange(GM, C) %>% pull(G)
	G1 = G1_linear %>% length
	G2_linear = counts_baseline %>% filter(!`house keeping`) %>% filter(level ==2) %>% distinct(G, GM, C) %>% arrange(GM, C) %>% pull(G)
	G2 = G2_linear %>% length
	G3_linear = counts_baseline %>% filter(!`house keeping`) %>% filter(level ==3) %>% distinct(G, GM, C) %>% arrange(GM, C) %>% pull(G)
	G3 = G3_linear %>% length


	# Observed mix counts
	y_linear_1 = y_source %>% filter(level ==1) %>% distinct(GM, Q, S, `read count`) %>% arrange(GM, Q) %>% pull(`read count`)
	y_linear_2 = y_source %>% filter(level ==2) %>% distinct(GM, Q, S, `read count`) %>% arrange(GM, Q) %>% pull(`read count`)
	y_linear_3 = y_source %>% filter(level ==3) %>% distinct(GM, Q, S, `read count`) %>% arrange(GM, Q) %>% pull(`read count`)

	# Observed mix samples indexes
	y_linear_S_1 = y_source %>% filter(level ==1) %>% distinct(GM, Q, S, `read count`) %>% arrange(GM, Q) %>% pull(S)
	y_linear_S_2 = y_source %>% filter(level ==2) %>% distinct(GM, Q, S, `read count`) %>% arrange(GM, Q) %>% pull(S)
	y_linear_S_3 = y_source %>% filter(level ==3) %>% distinct(GM, Q, S, `read count`) %>% arrange(GM, Q) %>% pull(S)

	# Observed mix samples indexes - get GM only for markers of each level. Exclude house keeping
	y_linear_GM_1 = y_source %>% filter(level ==1) %>% distinct(GM, Q, S, `read count`) %>% arrange(GM, Q) %>% pull(GM)
	y_linear_GM_2 = y_source %>% filter(level ==2) %>% distinct(GM, Q, S, `read count`) %>% arrange(GM, Q) %>% pull(GM)
	y_linear_GM_3 =y_source %>% filter(level ==3) %>% distinct(GM, Q, S, `read count`) %>% arrange(GM, Q) %>% pull(GM)

	# Lengths indexes
	Y_1 = y_linear_1 %>% length
	Y_2 = y_linear_2 %>% length
	Y_3 = y_linear_3 %>% length

	# Non centered
	lambda_mu_prior = c(6.2, 1)
	lambda_sigma_prior =  c( log(3.3) , 1)
	lambda_skew_prior =  c( -2.7, 1)
	sigma_intercept_prior = c( 1.9 , 0.1)
	lambda_log_scale = 	counts_baseline %>% filter(!query) %>% distinct(G, lambda) %>% arrange(G) %>% pull(lambda)

	# Horse shoe
	hs_df = 3             # If divergencies increase this
	par_ratio = 1/999     # number of expected non-zero divided by number of expected zeros
	hs_scale_slab = 2     # Expected value of the non-zeros
	df_global = 3;        # 1 == horseshoe, > 1 more like student-t
	df_slab = 25;         # stringency of the non-zero distribution. >> 1 if we are sure that non-zeros have exactly hs_scale_slab value, in a sense this increase the bimodality of the distribution

	# MODEL

	browser()

	fileConn<-file("~/.R/Makevars")
	writeLines(c( "CXX14FLAGS += -O3","CXX14FLAGS += -DSTAN_THREADS", "CXX14FLAGS += -pthread"), fileConn)
	close(fileConn)
	Sys.setenv("STAN_NUM_THREADS" = cores)
	ARMET_tc_model = stan_model("~/PhD/deconvolution/ARMET/src/stan_files/ARMET_tc.stan")

	Sys.time() %>% print
	fit =
		sampling(
			ARMET_tc_model, #stanmodels$ARMET_tc,
			chains=6, cores=6,
			iter=iterations, warmup=iterations-sampling_iterations,
			#include = F, pars=c("prop_a", "prop_b", "prop_c", "prop_d", "prop_e"),
			pars=c("prop_1", "prop_2", "prop_3", "exposure_rate", "error_ref_mix", "lambda_log", "sigma_inv_log", "sigma_intercept_dec"),
			#,
			init = function () list(	lambda_log = lambda_log_scale) # runif(G,  lambda_log_scale - 1, lambda_log_scale + 1)	)
			#save_warmup = FALSE,
			#pars = c("prop_1", "prop_2", "prop_3", "exposure_rate") #, "nb_sum") #,"mu_sum", "phi_sum"),
		) %>%
		{
			(.)  %>% summary() %$% summary %>% as_tibble(rownames="par") %>% arrange(Rhat %>% desc) %>% print
			(.)
		}
	Sys.time() %>% print

browser()
	ff = stan(file = "src/stan_files/horseshoe.stan")

	# # The reference fust 2 columns should look like this
	# counts_baseline %>% filter(!`house keeping`) %>% filter(level ==1) %>% distinct(`read count`, G, GM, C) %>% group_by(G, GM, C) %>% summarise(m = `read count` %>% median) %>% ungroup() %>% arrange(GM, C) %>% select(-G) %>% spread(GM, m) %>% select(2:3)
	#
	# fit %>% summary() %$% summary %>% as_tibble(rownames="par") %>%
	# 	separate(par, c(".variable", "C", "GM"), sep="[\\[,\\]]", extra="drop") %>%
	# 	mutate( C = C %>% as.integer, GM = GM %>% as.integer) %>% filter(grepl("mat_GM1", `.variable`))  %>%
	# 	left_join(
	# 		counts_baseline %>% filter(!`house keeping`) %>% filter(level ==1) %>% distinct(`read count`, G, GM, C) %>% group_by(G, GM, C) %>% summarise(m = `read count` %>% median) %>% ungroup() %>% arrange(GM, C)
	# 	) %>%
	# 	ggplot(aes(x=m+1, y=mean+1)) + geom_point() +scale_x_log10() + scale_y_log10()


	########################################
	# Parse results
	########################################

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
			tree %>% ToDataFrameTree("name", "C1", "C2", "C3") %>%
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
