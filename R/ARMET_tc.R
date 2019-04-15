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
						value_column = "read count normalised",
						log_transform = T
					) %>%

					# For immune cells repete twice
					{
						if((.) %>% distinct(`Cell type category`) %>% pull(1) == "immune_cell" )
							(.) %>%
							get_redundant_pair(
								distance_feature = "symbol",
								redundant_feature = "sample",
								value_column = "read count normalised",
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
	#reference = fit_df

# New mx
# mix_source =
# 	combn(ref_orig %>% distinct(C) %>% pull(1),m = 2) %>%
# 	t %>%
# 	as_tibble %>%
# 	mutate(run = 1:n()) %>%
# 	rowwise() %>%
# 	do({
# 		cc = (.)
# 		bind_rows(
# 			ref_orig %>% filter(C == cc %>% pull(1)) %>%
# 				sample_n(1) %>%
# 				distinct(sample, run),
# 			ref_orig %>% filter(C == cc %>% pull(2)) %>%
# 				sample_n(1) %>%
# 				distinct(sample, run)
# 		) %>%
# 			left_join(ref_orig) %>%
# 			distinct(`symbol`, `read count normalised`, `Cell type formatted`, run, sample)
# 	}) %>%
# 	ungroup()
#
# mix = mix_source %>%
# 	group_by(run) %>%
# 	distinct(`symbol`, `read count normalised`, `Cell type formatted`) %>%
# 	spread(`Cell type formatted`, `read count normalised`) %>%
# 	drop_na %>%
# 	mutate( `read count` = ( (`2` + `3`) / 2 ) %>% as.integer ) %>%
# 	mutate(sample = run) %>%
# 	unroup()

# old mix

mix_samples = c("ENCFF429MGN", "S00J8C11" )
mix_samples = c("counts.Fibroblast%20-%20Choroid%20Plexus%2c%20donor3.CNhs12620.11653-122E6", "C001FRB3" )
mix_samples = c("counts.Endothelial%20Cells%20-%20Microvascular%2c%20donor3.CNhs12024.11414-118F1","counts.CD4%2b%20T%20Cells%2c%20donor1.CNhs10853.11225-116C1" )

	mix =
		ref %>%
		inner_join( (.) %>% distinct(sample) %>% filter(sample %in% mix_samples)) %>%
		distinct(`symbol`, `read count normalised`, `Cell type formatted`) %>%
		spread(`Cell type formatted`, `read count normalised`) %>%
		drop_na %>%
		mutate( `read count` = ( (endothelial + t_CD4) / 2 ) %>% as.integer ) %>%
		mutate(sample = "1") %>%
		select(-c(2:3)) %>%
		spread(`symbol`, `read count`)

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

	# mix = reference %>%
	# 	inner_join( (.) %>% distinct(sample) %>% slice(c(200))) %>%
	# 	distinct(sample, `symbol`, `read count`) %>%
	# 	spread(`symbol`, `read count`) %>%
	# 	mutate(sample = 1:n() %>% as.character)

	house_keeping =
		ref %>% filter(`Cell type category` == "house_keeping" ) %>%
		distinct(`symbol`) %>%
		pull(1) %>%
		head(n=200)

	load("/wehisan/home/allstaff/m/mangiola.s/PhD/deconvolution/ARMET/data/tree_Yaml.rda")
	source("R/ARMET_parseTree.R")
	my_tree = format_tree( get_node_from_name(tree, "TME"), mix, "")
	markers = get_node_label_level_specfic(
		get_node_from_name(my_tree,  "TME"),
		label = "markers",
		start_level = 1,
		stop_level = 1
	)
	# markers =  read_csv("docs/markers.csv") %>% pull(symbol)

	reference =
		ref %>%

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
		) %>%

		# Add cell type indeces
		left_join(
			(.) %>%
				filter(!`Cell type category` %in% c("query", "house_keeping")  ) %>%
				distinct(`Cell type category`) %>%
				mutate(C = 1:n())
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

	y =
		df %>%
		filter(`Cell type category` == "query") %>%
		select(S, Q, `symbol`, `read count`) %>%
		left_join(
			counts_baseline %>%
				distinct(`symbol`, G, `Cell type category`)
		) %>%
		spread(`Cell type category`, G) %>%
		select(-`symbol`) %>%
		select(`read count`, S, Q, everything()) %>%
		as_matrix

	I = y %>% nrow
	Q = df %>% filter(`Cell type category` == "query") %>% distinct(Q) %>% nrow



	fileConn<-file("~/.R/Makevars")
	writeLines(c( "CXX14FLAGS += -O3","CXX14FLAGS += -DSTAN_THREADS", "CXX14FLAGS += -pthread"), fileConn)
	close(fileConn)
	Sys.setenv("STAN_NUM_THREADS" = cores)
	ARMET_tc = stan_model("src/stan_files/ARMET_tc.stan")


	initializer <- function() list(
		lambda_log =
			counts_baseline %>%
			# Ths is bcause mix lacks lambda info and produces NA in the df
			filter(!(`Cell type category` == "house_keeping" & lambda %>% is.na)) %>%
			distinct(G, lambda) %>% pull(lambda),
		sigma_raw =
			counts_baseline %>%
			# Ths is bcause mix lacks lambda info and produces NA in the df
			filter(!(`Cell type category` == "house_keeping" & sigma_raw %>% is.na)) %>%
			distinct(G, sigma_raw) %>% pull(sigma_raw)
	)
	inits <- foreach(chain = 1:3) %do% {initializer()}


	Sys.time()
	fit =
		sampling(
			ARMET_tc, #stanmodels$ARMET_tc,
			chains=3, cores=3,
			iter=300, warmup=200,init=inits
			# ,
			# save_warmup = FALSE,
			# pars = c(
			# 	"lambda", "sigma_raw", "sigma_intercept", "sigma_slope",
			# 	"sigma_sigma", "lambda_mu", "lambda_sigma", "lambda_skew"
			# )
		)
	Sys.time()

	# Plot differences in lambda
	 (fit %>%
		gather_draws(lambda_log[G]) %>%
		median_qi() %>%
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




	writeLines("ARMET: building output")

	# # Create tree with hypothesis testing
	# library(data.tree)
	# fileName=("~/PhD/deconvolution/ARMET_dev/ARMET_TME_tree.yaml")
	# ya=readChar(fileName, file.info(fileName)$size)
	# osList <- yaml::yaml.load(ya)
	# treeYaml <- as.Node(osList)
	# save(treeYaml, file="data/tree_Yaml.rda")

	treeYaml.stat =
		switch(
			(!is.null(cov_to_test)) + 1,
			NULL,
			get_tree_hypoth_test(treeYaml, my_tree, cov_to_test)
		)

	if(save_report) save(treeYaml.stat, file=sprintf("%s/tree_pvalues.RData", output_dir))

	# Return
	list(

		# Matrix of proportions
		proportions =	get_last_existing_leaves_with_annotation( my_tree ) %>%
			select(-relative_proportion),

		# What mixture was used by the model after normalization
		mix = mix %>%
			filter(gene %in% get_genes( my_tree )) %>%
			droplevels(),

		# Return the statistics
		stats = treeYaml.stat,

		# Return the annotated tree
		tree = my_tree,

		# Return the input itself
		input = input
	)

}
