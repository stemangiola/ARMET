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
	source("https://gist.githubusercontent.com/stemangiola/90a528038b8c52b21f9cfa6bb186d583/raw/9f361c5ecd5e15473e831945971cf290ea7377f1/transcription_tool_kit.R")
	library(tidyverse)
	library(foreach)
	library(rstan)
	cores = 14
	input = c(as.list(environment()))
	shards = 56

	# Global properties
	sigma_intercept = 1.3663737
	sigma_slope = -0.3443049
	sigma_sigma = 1.1755951
	lambda_mu_mu = 7
	lambda_skew = -8.7458016
	lambda_sigma = 8.7234045

	format_for_MPI = function(df){
		df %>%

		left_join(
			(.) %>%
				distinct(symbol) %>%
				mutate( idx_MPI = head( rep(1:shards, (.) %>% nrow %>% `/` (shards) %>% ceiling ), n=(.) %>% nrow) )
		) %>%
		arrange(idx_MPI, symbol) %>%

		# Decide start - end location
		group_by(idx_MPI) %>%
		do(
			(.) %>%
				left_join(
					(.) %>%
						distinct(idx_MPI, sample, symbol) %>%
						arrange(idx_MPI, symbol) %>%
						count(idx_MPI, symbol) %>%
						mutate(end = cumsum(n)) %>%
						mutate(start = c(1, .$end %>% rev() %>% `[` (-1) %>% rev %>% `+` (1)))
				)
		) %>%
		ungroup() %>%

		# Add counts MPI rows indexes
		group_by(idx_MPI) %>%
		mutate(`read count MPI row` = 1:n()) %>%
		ungroup %>%

		# Add symbol MPI rows indexes
		left_join(
			(.) %>%
				group_by(idx_MPI) %>%
				distinct(symbol) %>%
				mutate(`symbol MPI row` = 1:n()) %>%
				ungroup
		) %>%

		# Add gene idx
		left_join(
			(.) %>%
				distinct(symbol, idx_MPI, `symbol MPI row`) %>%
				arrange(idx_MPI, `symbol MPI row`) %>%
				mutate(G = 1:n())
		)

	}

	######################################
	# Variable should be already set
	######################################

	load("data/reference.RData")

	house_keeping =
		reference %>% filter(`Cell type category` == "house_keeping" ) %>%
		distinct(`symbol original`) %>%
		pull(1) %>%
		head(n=200)

	reference =
		reference %>%

		filter(
			#`Cell type category` == "house_keeping" |
			`symbol original` %in% 	house_keeping |
			`symbol original` %in%  ( read_csv("docs/markers.csv") %>% pull(symbol))
		) %>%
		select( -start, -end, -n, -contains("idx")) %>%
		mutate(`read count` = `read count` %>% as.integer)
		#inner_join( (.) %>% distinct(sample) %>% head(n=100))

	# reference %>%
	# 	# Reduce numbers by redundancy
	# 	add_MDS_components(
	# 		replicates_column = "symbol original",
	# 		cluster_by_column = "sample",
	# 		value_column = "read count normalised"
	# 	)


	mix = reference %>%
		inner_join( (.) %>% distinct(sample) %>% head(n=1)) %>%
		distinct(sample, `symbol original`, `read count`) %>%
		spread(`symbol original`, `read count`) %>%
		mutate(sample = 1:n() %>% as.character)

	######################################
	######################################

	# Merge data sets
	df =
		bind_rows(
			reference %>%
			filter(`symbol original` %in% ( mix %>% colnames )),
			mix %>%
				gather(`symbol original`, `read count`, -sample) %>%
				mutate(`Cell type category` = ifelse(
					`symbol original` %in% ( reference %>% filter(`Cell type category` == "house_keeping" ) %>% distinct(`symbol original`) %>% pull(1) 	),
					"house_keeping",
					"query"
				)) %>%
				unite(symbol, c("Cell type category", "symbol original"), remove = F)
		)	%>%

		# Add symbol indeces
		left_join(
			(.) %>%
				filter(!`Cell type category` %in% c("house_keeping")  ) %>%
				distinct(`symbol original`) %>%
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
		filter(`Cell type category` != "query") %>%
		format_for_MPI

	N = counts_baseline %>% distinct(sample, symbol, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% max
	M = counts_baseline %>% distinct(start, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% max
	G = counts_baseline %>% distinct(symbol) %>% nrow()
	S = counts_baseline %>% distinct(sample) %>% nrow()
	G_per_shard = counts_baseline %>% distinct(symbol, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% as.array
	n_shards = min(shards, counts_baseline %>% distinct(idx_MPI) %>% nrow)
	G_per_shard_idx = c(0, counts_baseline %>% distinct(symbol, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% cumsum)

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
		select(S, Q, `symbol original`, `read count`) %>%
		left_join(
			counts_baseline %>%
				distinct(`symbol original`, G, C)
		) %>%
		spread(C, G) %>%
		select(-`symbol original`) %>%
		select(`read count`, S, Q, everything()) %>%
		as_matrix

	I = y %>% nrow
	Q = df %>% filter(`Cell type category` == "query") %>% distinct(Q) %>% nrow



	fileConn<-file("~/.R/Makevars")
	writeLines(c("CXX14 = g++ ", "CXX14FLAGS += -O3","CXX14FLAGS += -DSTAN_THREADS", "CXX14FLAGS += -pthread"), fileConn)
	close(fileConn)
	Sys.setenv("STAN_NUM_THREADS" = cores)
	ARMET_tc = stan_model("src/stan_files/ARMET_tc.stan")

	Sys.time()
	fit =
		sampling(
			ARMET_tc, #stanmodels$ARMET_tc,
			chains=3, cores=3,
			iter=300, warmup=200
			# ,
			# save_warmup = FALSE,
			# pars = c(
			# 	"lambda", "sigma_raw", "sigma_intercept", "sigma_slope",
			# 	"sigma_sigma", "lambda_mu", "lambda_sigma", "lambda_skew"
			# )
		)
	Sys.time()


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
