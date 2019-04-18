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
	ct_to_omit =                        c("t_CD4_naive", "adipocyte"),
	verbose =                           F,
	omit_regression =                   F,
	save_fit =                          F,
	seed =                              NULL,
	cores = 14
){

	source("https://gist.githubusercontent.com/stemangiola/dd3573be22492fc03856cd2c53a755a9/raw/e4ec6a2348efc2f62b88f10b12e70f4c6273a10a/tidy_extensions.R")
	source("https://gist.githubusercontent.com/stemangiola/90a528038b8c52b21f9cfa6bb186d583/raw/5d95739c87a26568f96f8dc7d29528473fbe82c2/transcription_tool_kit.R")
	library(tidyverse)
	library(foreach)
	library(rstan)
	cores = 14
	input = c(as.list(environment()))
	shards = 56

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

	# Global properties
	sigma_intercept = 1.3420415
	sigma_slope = -0.3386389
	sigma_sigma = 1.1720851
	lambda_mu_mu = 5.612671
	lambda_sigma = 7.131593

	format_for_MPI = function(df){
		df %>%

			left_join(
				(.) %>%
					distinct(ct_symbol) %>%
					mutate( idx_MPI = head( rep(1:shards, (.) %>% nrow %>% `/` (shards) %>% ceiling ), n=(.) %>% nrow) )
			) %>%
			arrange(idx_MPI, ct_symbol) %>%

			# Decide start - end location
			group_by(idx_MPI) %>%
			do(
				(.) %>%
					left_join(
						(.) %>%
							distinct(idx_MPI, sample, ct_symbol) %>%
							arrange(idx_MPI, ct_symbol) %>%
							count(idx_MPI, ct_symbol) %>%
							mutate(end = cumsum(n)) %>%
							mutate(start = c(1, .$end %>% rev() %>% `[` (-1) %>% rev %>% `+` (1)))
					)
			) %>%
			ungroup() %>%

			# Add counts MPI rows indexes
			group_by(idx_MPI) %>%
			mutate(`read count MPI row` = 1:n()) %>%
			ungroup %>%

			# Add ct_symbol MPI rows indexes
			left_join(
				(.) %>%
					group_by(idx_MPI) %>%
					distinct(ct_symbol) %>%
					mutate(`symbol MPI row` = 1:n()) %>%
					ungroup
			) %>%

			# Add gene idx
			left_join(
				(.) %>%
					distinct(ct_symbol, idx_MPI, `symbol MPI row`) %>%
					arrange(idx_MPI, `symbol MPI row`) %>%
					mutate(G = 1:n())
			)

	}

	decrease_replicates = function(my_df){

		my_df %>%

			# Correct house keeping formatting
			mutate(`house keeping` = `Cell type category` == "house_keeping") %>%
			rename(`Cell type category old` = `Cell type category`) %>%
			left_join(
				(.) %>%
					distinct(sample, `Cell type category old`) %>%
					filter(`Cell type category old` != "house_keeping") %>%
					rename(`Cell type category` = `Cell type category old`)
			) %>%

			# Detect redundant
			multidplyr::partition(`Cell type category`, `Data base`) %>%
			#group_by(`Cell type category`, `Data base`) %>%
			do({
				`%>%` = magrittr::`%>%`
				library(tidyverse)
				library(magrittr)
				library(foreach)
				source("https://gist.githubusercontent.com/stemangiola/dd3573be22492fc03856cd2c53a755a9/raw/e4ec6a2348efc2f62b88f10b12e70f4c6273a10a/tidy_extensions.R")
				source("https://gist.githubusercontent.com/stemangiola/90a528038b8c52b21f9cfa6bb186d583/raw/4ac3a7d74d4f7971bd8dd8eb470279eb7f42bbf8/transcription_tool_kit.R")

				(.) %>%
					get_redundant_pair(
						distance_feature = "symbol",
						redundant_feature = "sample",
						value_column = "read count normalised bayes",
						log_transform = T
					) %>%

					# For immune cells repete twice
					{
						if((.) %>% distinct(`Cell type category`) %>% pull(1) == "immune_cell" )
							(.) %>%
							get_redundant_pair(
								distance_feature = "symbol",
								redundant_feature = "sample",
								value_column = "read count normalised bayes",
								log_transform = T
							)
						else
							(.)
					}

			}) %>%

			collect %>% ungroup %>%
			# Re put house keeping
			mutate(`Cell type category` = ifelse(`house keeping`, "house_keeping", `Cell type category`))

	}

	######################################
	# Variable should be already set
	######################################

	load("docs/ref.RData")
	sample_blacklist = c("666CRI", "972UYG", "344KCP", "555QVG", "370KKZ", "511TST", "13816.11933", "13819.11936", "13817.11934", "13818.11935", "096DQV", "711SNV")

# New mx
set.seed(123)
mix_source =
	combn(ref %>% filter(`Cell type category` != "house_keeping") %>% distinct(`Cell type category`) %>% pull(1),m = 2) %>%
	t %>%
	as_tibble %>%
	mutate(run = 1:n()) %>%
	group_by(run) %>%
	do({
		cc = (.)
		bind_rows(
			ref %>%
				filter( grepl(sample_blacklist %>% paste(collapse="|"), sample)) %>%
				filter(`Cell type category` == (cc %>% pull(1))) %>%
				sample_n(1) %>%
				distinct(sample, run),
			ref %>%
				filter( grepl(sample_blacklist %>% paste(collapse="|"), sample)) %>%
				filter(`Cell type category` == (cc %>% pull(2))) %>%
				sample_n(1) %>%
				distinct(sample, run)
		) %>%
			left_join(ref) %>%
			distinct(`symbol`, `read count normalised bayes`, `Cell type formatted`, run, sample)
	}) %>%
	ungroup()

mix = mix_source %>%
	group_by(run) %>%
	do(
		(.) %>%
			distinct(`symbol`, `read count normalised bayes`, `Cell type formatted`, run) %>%
			spread(`Cell type formatted`, `read count normalised bayes`) %>%
			drop_na %>%
			mutate(combination = names((.))[3:4] %>% paste(collapse=" ")) %>%
			setNames(c("symbol", "run", "1", "2", "combination")) %>%
			mutate( `read count` = ( (`1` + `2`) / 2 ) %>% as.integer ) %>%
			unite(sample, c("run", "combination"), remove = F)
	) %>%
	ungroup() %>%
	select(sample, symbol, `read count`) %>%
	spread(`sample`, `read count`) %>%
	drop_na %>%
	gather(sample, `read count`, -symbol) %>%
	spread(`symbol`, `read count`)

# load("/wehisan/home/allstaff/m/mangiola.s/PhD/deconvolution/ARMET/data/tree_Yaml.rda")
# source("R/ARMET_parseTree.R")
# my_tree = format_tree( get_node_from_name(tree, "TME"), mix, "")
# markers = get_node_label_level_specfic(
# 	get_node_from_name(my_tree,  "TME"),
# 	label = "markers",
# 	start_level = 1,
# 	stop_level = 1
# )
markers =  read_csv("docs/markers.csv") %>% pull(symbol) %>% unique


# Trouble shoot
plot_df =
	ref %>%
	filter(filt_for_calc) %>%
	decrease_replicates %>%
	filter(	`symbol` %in%  markers ) %>%
	bind_rows(
		mix %>%
			gather(symbol, `read count`, -sample) %>%
			mutate(`Cell type category` = "mix") %>%
			mutate(`read count normalised bayes` = `read count`)
	) %>%
	add_MDS_components(
		replicates_column = "symbol",
		cluster_by_column = "sample",
		value_column = "read count normalised bayes",
		log_transform = T
	)

plot_df %>%
	filter(`Cell type category` != "house_keeping") %>%
	distinct(`PC 1`, `PC 2`, `Cell type category`, `Data base`, sample) %>%
	{
		(.) %>% ggplot(aes(
		x=`PC 1`, y = `PC 2`,
		color=`Cell type category`,
		label=sample, ds=`Data base`
		)) +
		geom_point() + my_theme
	} %>%
	plotly::ggplotly()



	# ar = ARMET_tc(
	# 	ref_orig %>%
	# 		inner_join( (.) %>% distinct(sample) %>% filter(sample %in% mix_samples)) %>%
	# 		distinct(`symbol`, `read count normalised`, `Cell type formatted`) %>%
	# 		spread(`Cell type formatted`, `read count normalised`) %>%
	# 		drop_na %>%
	# 		mutate( `read count` = ( (macrophage_M2 + endothelial) / 2 ) %>% as.integer ) %>%
	# 		mutate(sample = "1") %>%
	# 		select(-c(2:3)) %>% rename(gene = `symbol`) %>% spread(sample, `read count`), do_debug=F
	# )

	house_keeping =
		ref %>% filter(`Cell type category` == "house_keeping" ) %>%
		distinct(`symbol`) %>%
		pull(1) %>%
		head(n=200)



	reference =
		ref %>%
		filter( grepl(sample_blacklist %>% paste(collapse="|"), sample)) %>%

		# Decrese number of samples
		decrease_replicates %>%

		# Get only markes and house keeping
		filter(
			#`Cell type category` == "house_keeping" |
			`symbol` %in% 	house_keeping |
				`symbol` %in%  markers # %>% sample %>% head(n=100))
		) %>%
		select(  -contains("idx"), -G) %>%
		mutate(`read count` = `read count` %>% as.integer)

	######################################
	######################################

	# Merge data sets
	df =
		bind_rows(
			reference %>%
				filter(`symbol` %in% ( mix %>% colnames )),
			mix %>%
				gather(`symbol`, `read count`, -sample) %>%
				inner_join(reference %>% distinct(symbol) ) %>%
				mutate(`Cell type category` = ifelse(
					`symbol` %in% ( reference %>% filter(`Cell type category` == "house_keeping" ) %>% distinct(`symbol`) %>% pull(1) 	),
					"house_keeping",
					"query"
				))
		)	%>%

		# Add symbol indeces
		left_join(
			(.) %>%
				filter(!`Cell type category` %in% c("house_keeping")  ) %>%
				distinct(`symbol`) %>%
				mutate(M = 1:n())
		) %>%

		# Add sample indeces
		mutate(S = sample %>% as.factor %>% as.integer) %>%
		left_join(
			(.) %>%
				filter(`Cell type category` == "query"  ) %>%
				distinct(`sample`) %>%
				mutate(Q = 1:n())
		)

	# For  reference MPI inference
	counts_baseline =
		df %>%
		filter(`Cell type category` != "query")  %>%
		unite(ct_symbol, c("Cell type category", "symbol"), remove = F) %>%

		format_for_MPI

	N = counts_baseline %>% distinct(sample, ct_symbol, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% max
	M = counts_baseline %>% distinct(start, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% max
	G = counts_baseline %>% distinct(ct_symbol) %>% nrow()
	S = counts_baseline %>% distinct(sample) %>% nrow()
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

	# For deconvolution
	ct_in_levels = c(4) %>% as.array

	y_source =
		df %>%
		filter(`Cell type category` == "query") %>%
		select(S, Q, `symbol`, `read count`) %>%
		left_join(
			counts_baseline %>%
				distinct(`symbol`, G, `Cell type category`)
		)
	y =
		y_source %>%
		spread(`Cell type category`, G) %>%
		select(-`symbol`) %>%
		select(`read count`, S, Q, everything()) %>%
		as_matrix

	I = y %>% nrow
	Q = df %>% filter(`Cell type category` == "query") %>% distinct(Q) %>% nrow

	# Pass previously infered parameters
	do_infer = 0
	lambda_log_data =
		counts_baseline %>%
		# Ths is bcause mix lacks lambda info and produces NA in the df
		filter(!(`Cell type category` == "house_keeping" & lambda %>% is.na)) %>%
		distinct(G, lambda) %>% pull(lambda)

	sigma_raw_data =
		counts_baseline %>%
		# Ths is bcause mix lacks lambda info and produces NA in the df
		filter(!(`Cell type category` == "house_keeping" & sigma_raw %>% is.na)) %>%
		distinct(G, sigma_raw) %>% pull(sigma_raw)


	fileConn<-file("~/.R/Makevars")
	writeLines(c( "CXX14FLAGS += -O3","CXX14FLAGS += -DSTAN_THREADS", "CXX14FLAGS += -pthread"), fileConn)
	close(fileConn)
	Sys.setenv("STAN_NUM_THREADS" = cores)
	ARMET_tc = stan_model("src/stan_files/ARMET_tc.stan")

	Sys.time()
	fit =
		sampling(
			ARMET_tc, #stanmodels$ARMET_tc,
			chains=3, cores=3,
			iter=300, warmup=200
			#,
			#init=foreach(chain = 1:3) %do% {(function() list( lambda_log = lambda_log_data, sigma_raw = sigma_raw_data	))()}
			# ,
			# save_warmup = FALSE,
			# pars = c(
			# 	"lambda", "sigma_raw", "sigma_intercept", "sigma_slope",
			# 	"sigma_sigma", "lambda_mu", "lambda_sigma", "lambda_skew"
			# )
		)
	Sys.time()

	# Produce results
	prop = fit %>%
		tidybayes::gather_draws(prop[Q, C]) %>%
		tidybayes::median_qi() %>%
		left_join(
			y_source %>%
				distinct(`Cell type category`) %>%
				arrange(`Cell type category`) %>%
				mutate(C = 1:n())
		) %>%
		left_join(
			df %>%
		 	filter(`Cell type category` == "query") %>%
		 	select(Q, sample)
		)

if(0){

	# Plot differences in lambda
	 (fit %>%
		tidybayes::gather_draws(lambda_log[G]) %>%
		tidybayes::median_qi() %>%
		left_join(
			counts_baseline %>%
				distinct(`symbol`, G, `Cell type category`)
		) %>%
		left_join(
			reference %>%
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


	# Return
	list(

		# Matrix of proportions
		proportions =	prop,

		# Return the input itself
		input = input
	)

}
