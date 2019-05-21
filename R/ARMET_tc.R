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
get_MPI_deconv = function(y_source, shards, my_level){
	y_MPI_source =
		y_source %>%
		distinct(level, Q, S, symbol, G, `Cell type category`, `read count`) %>%
		filter(level==my_level)  %>%
		arrange(Q, symbol) %>%
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

get_plot_predicted_real = function(fit_parsed, y_source){
	fit_parsed %>%
		tidybayes::median_qi() %>%
		bind_cols(
			y_source %>% distinct(level, Q, S, symbol, `read count`) %>% arrange(level, Q, symbol)
		) %>%
		select(-mu_sum.lower, -mu_sum.upper, -phi_sum.lower, -phi_sum.upper) %>%
		group_by(I) %>%
		do(
			bind_cols(
				(.),
				rnbinom(10000, mu = (.) %>% pull(mu_sum), size = (.) %>% pull(phi_sum)) %>% quantile(probs = c(0.1, 0.9)) %>% enframe %>% spread(name, value)
			)
		) %>%
		ungroup() %>%
		rowwise %>%
		mutate(inside = between(`read count`, `10%`, `90%`)) %>%
		ungroup %>%
		ggplot(aes(x=`read count` + 1, y=mu_sum + 1, color=inside)) +
		geom_point(alpha=0.5)  +
		geom_abline(slope = 1, intercept = 0) +
		geom_errorbar(aes(ymin=`10%`, ymax=`90%`), alpha=0.5) +
		facet_wrap(~Q) +
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
									slice(1: n_markers)
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
	iterations = 300
){

	# full_bayesian = 0
	# ct_to_omit =                        c("t_CD4_naive", "adipocyte")
	# verbose =                           F
	# omit_regression =                   F
	# save_fit =                          F
	# seed =                              NULL
	# cores = 14

	source("https://gist.githubusercontent.com/stemangiola/dd3573be22492fc03856cd2c53a755a9/raw/e4ec6a2348efc2f62b88f10b12e70f4c6273a10a/tidy_extensions.R")
	source("https://gist.githubusercontent.com/stemangiola/90a528038b8c52b21f9cfa6bb186d583/raw/4a5798857362d946bd3029188b1cc9eb9b625456/transcription_tool_kit.R")
	library(tidyverse)
	library(foreach)
	library(rstan)

	input = c(as.list(environment()))
	shards = cores * 4

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

	# Print overlap descriptive stats
	#get_overlap_descriptive_stats(mix %>% slice(1) %>% gather(`symbol`, `read count`, -sample), reference)

	# Merge data sets


	reference_filtered = filter_reference(reference, mix)

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
			(.) %>% filter(`house keeping` & !`query`) %>% distinct(symbol, level) %>% group_by(symbol) %>% arrange(level) %>% slice(2) %>% ungroup()
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
		)

	G = df %>% filter(!`query`) %>% distinct(G) %>% nrow()

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

	ct_in_levels = c(4,7)

	y_source =
		df %>%
		filter(`query` & !`house keeping`) %>%
		select(S, Q, `symbol`, `read count`) %>%
		left_join(	df %>% filter(!query) %>% distinct(`symbol`, G, `Cell type category`, level, lambda, sigma_raw) ) %>%
		arrange(level, `Cell type category`, Q, symbol) %>%
		mutate(`Cell type category` = factor(`Cell type category`, unique(`Cell type category`)))

	idx_1_source = y_source %>% filter(level == 1) %>% distinct(symbol, G, `Cell type category`) %>% arrange(`Cell type category`, symbol)
	idx_2_source = y_source %>% filter(level == 2) %>% distinct(symbol, G, `Cell type category`) %>% arrange(`Cell type category`, symbol)
	idx_1 = idx_1_source %>% pull(G)
	idx_2 = idx_2_source %>% pull(G)
	I1 = idx_1 %>% length
	I2 = idx_2 %>% length
	I1_dim = c(idx_1_source %>% distinct(symbol) %>% nrow, idx_1_source %>% distinct(`Cell type category`) %>% nrow)
	I2_dim = c(idx_2_source %>% distinct(symbol) %>% nrow, idx_2_source %>% distinct(`Cell type category`) %>% nrow)

	# Data MPI for deconvolution level 1
	y_MPI_lv1 = y_source %>% get_MPI_deconv(shards, 1)
	y_MPI_source_lv1 = y_MPI_lv1 %$% y_MPI_source
	y_MPI_symbol_per_shard_lv1 = y_MPI_lv1 %$% y_MPI_symbol_per_shard
	y_MPI_G_per_shard_lv1 = y_MPI_lv1 %$% y_MPI_G_per_shard
	y_MPI_idx_lv1 = y_MPI_lv1 %$% y_MPI_idx
	y_MPI_N_per_shard_lv1 = y_MPI_lv1 %$% y_MPI_N_per_shard
	y_MPI_count_lv1 = y_MPI_lv1 %$% y_MPI_count

	# Data MPI for deconvolution level 2
	y_MPI_lv2 = y_source %>% get_MPI_deconv(shards, 2)
	y_MPI_source_lv2 = y_MPI_lv2 %$% y_MPI_source
	y_MPI_symbol_per_shard_lv2 = y_MPI_lv2 %$% y_MPI_symbol_per_shard
	y_MPI_G_per_shard_lv2 = y_MPI_lv2 %$% y_MPI_G_per_shard
	y_MPI_idx_lv2 = y_MPI_lv2 %$% y_MPI_idx
	y_MPI_N_per_shard_lv2 = y_MPI_lv2 %$% y_MPI_N_per_shard
	y_MPI_count_lv2 = y_MPI_lv2 %$% y_MPI_count

	Q = df %>% filter(`query`) %>% distinct(Q) %>% nrow
	idx_ct_root = c(1:4)
	idx_ct_immune = c(1:3, 5:11)
	y_idx_ct_root = idx_ct_root + 3
	y_idx_ct_immune = idx_ct_immune + 3

	#######################################
	# Merge all MPI
	#######################################

	# Reference
	counts_package =
		rep(c(M, N, S), shards) %>% matrix(nrow = shards, byrow = T) %>%
		cbind(G_per_shard) %>%
		cbind(symbol_end) %>%
		cbind(sample_idx) %>%
		cbind(counts)

	# level 1
	lev1_package =
		rep(c(ct_in_levels[1], Q, S), shards) %>% matrix(nrow = shards, byrow = T) %>%
		cbind(y_MPI_symbol_per_shard_lv1) %>%
		cbind(y_MPI_G_per_shard_lv1) %>%
		cbind(y_MPI_N_per_shard_lv1) %>%
		cbind(y_MPI_count_lv1)

	# level 1
	lev2_package =
		rep(c(sum(ct_in_levels[1:2]) - 1, Q, S), shards) %>% matrix(nrow = shards, byrow = T) %>%
		cbind(y_MPI_symbol_per_shard_lv2) %>%
		cbind(y_MPI_G_per_shard_lv2) %>%
		cbind(y_MPI_N_per_shard_lv2) %>%
		cbind(y_MPI_count_lv2)

	# Integrate everything
	data_package =

		# Size 3 data data sets
		rep(c(ncol(counts_package), ncol(lev1_package), ncol(lev2_package)), shards) %>%
		matrix(nrow = shards, byrow = T) %>%

		# Size parameter datasets
		cbind(
			rep(c(
				(2*M + S),
				(max(y_MPI_G_per_shard_lv1) * 2 + Q + Q + (Q * ct_in_levels[1])),
				(max(y_MPI_G_per_shard_lv2) * 2 + Q + Q + (Q * (sum(ct_in_levels) - 1)))
			), shards) %>%
			matrix(nrow = shards, byrow = T)
		) %>%

		# Data sets
		cbind(counts_package) %>%
		cbind(lev1_package) %>%
		cbind(lev2_package)

	########################################
	########################################

	# Old daa structure

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

	# Testing
	# exposure_rate = df %>% distinct(S) %>% nrow %>% seq(-1, 1, length.out = .);
	# set.seed(143)
	# prop_1 = gtools::rdirichlet(Q, c(1,1,1,1))

	# browser()
	# load("temp_fit.RData")
	# exposure_rate = c(1.5, 1.6, 1.9, 1.8, 1.5, 1.9, 2.0, 1.7, 1.6, 1.5)

	fileConn<-file("~/.R/Makevars")
	writeLines(c( "CXX14FLAGS += -O3","CXX14FLAGS += -DSTAN_THREADS", "CXX14FLAGS += -pthread"), fileConn)
	close(fileConn)
	Sys.setenv("STAN_NUM_THREADS" = cores)
	ARMET_tc = stan_model("src/stan_files/ARMET_tc.stan")



	Sys.time() %>% print
	fit =
		sampling(
			ARMET_tc, #stanmodels$ARMET_tc,
			chains=3, cores=3,
			iter=iterations, warmup=iterations-100,
			pars = c("prop_1", "prop_2", "exposure_rate", "sigma_raw_global") #,"mu_sum", "phi_sum")
		)
	Sys.time() %>% print


	#browser()

	# Produce results
	prop =
		fit %>%
		tidybayes::gather_draws(prop_1[Q, C], prop_2[Q, C]) %>%
		filter(.variable %in% c("prop_1", "prop_2")) %>%

		# If not converged choose the majority chains
		mutate(	converged = diptest::dip.test(`.value`) %$%	`p.value` > 0.05) %>%

		do({
			if((.) %>% distinct(converged) %>% pull(1)) (.)
			else (.) %>% choose_chains_majority_roule
		}) %>%

		# Summarise
		{
			left_join(
				(.) %>% select(-converged) %>% tidybayes::median_qi(),
				(.) %>% distinct(converged)
			)
		} %>%
		ungroup() %>%

		# Parse
		separate(.variable, c(".variable", "level"), convert = T) %>%
		left_join(

			y_source %>%
				distinct(`Cell type category`) %>%
				arrange(`Cell type category`) %>%
				{
					ys = (.)
					ys %>%
						slice(!!idx_ct_root) %>%
						mutate(C = 1:n(), level=1) %>%
						bind_rows(
							ys %>%
								slice(!!idx_ct_immune) %>%
								mutate(C = 1:n(), level=2)
						)
				}
		) %>%
		left_join(
			df %>%
				filter(`query`) %>%
				distinct(Q, sample)
		)


	#get_plot_predicted_real( fit %>% tidybayes::spread_draws(mu_sum[I], phi_sum[I]) %>% filter(.chain==2), y_source ) + my_theme


	# Return
	list(

		# Matrix of proportions
		proportions =	prop,

		# Return the input itself
		input = input,

		# Return the fitted object
		fit = fit
	)

}
