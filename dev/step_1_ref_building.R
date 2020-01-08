# Infer NB from Cell type category

# Initialise
#setwd("/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/RNAseq-noise-model")
set.seed(13254)

# Import libraries
library(tidyverse)
library(magrittr)
library(rstan)
library(tidybayes)
library(here)
library(foreach)
# library(parallel)
# library(doParallel)
library(multidplyr)
library(data.tree)
library(ttBulk)

# library(future)
# plan(multiprocess)

# rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())

# Get level
args <- commandArgs(TRUE)
my_level = args[1] %>% as.integer

# Registering parallel framework
system("nproc", intern = TRUE) %>%
	as.integer %>%
	sprintf("Working on %s cores", .) %>%
	print

n_cores = system("nproc", intern = TRUE) %>% as.integer %>% sum(-4)


# cl = n_cores %>% makeCluster( manual=FALSE, outfile='log.txt')
# clusterEvalQ(cl, library(tidyverse))
# clusterEvalQ(cl, library(magrittr))
# registerDoParallel(cl)
#set_default_cluster(cl)

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

# Tools
source("../../../PostDoc/ppcSeq/R/do_parallel.R")
source("https://gist.githubusercontent.com/stemangiola/9d2ba5d599b7ac80404c753cdee04a01/raw/26e5b48fde0cd4f5b0fd7cbf2fde6081a5f63e7f/tidy_data_tree.R")

ifelse_pipe = function(.x, .p, .f1, .f2 = NULL) {
	switch(.p %>% `!` %>% sum(1),
				 as_mapper(.f1)(.x),
				 if (.f2 %>% is.null %>% `!`)
				 	as_mapper(.f2)(.x)
				 else
				 	.x)

}


#options(error = quote({dump.frames(to.file=TRUE); q()}))

# Setup table of name conversion
tree =
	yaml:: yaml.load_file("/wehisan/home/allstaff/m/mangiola.s/PhD/deconvolution/ARMET/data/tree.yaml") %>%
	data.tree::as.Node()


# MPI
do_inference = function(counts, my_run){

	shards = n_cores %>% divide_by(3) %>% floor %>% multiply_by(4)

	counts = 	counts %>%

		# select run and HKG
		filter(grepl("house_keeping", symbol) | run == my_run)
		#filter(run == my_run)


	res =
		counts %>%
		mutate(do_check = `Cell type category` == "house keeping") %>%
		mutate(PValue = ifelse(do_check, 0, 1)) %>%
		mutate(`count` = `count` %>% as.integer) %>%
		inner_join((.) %>% distinct(symbol) %>% sample_frac(0.1)) %>%

		ppcSeq::ppc_seq(
			significance_column = PValue,
			do_check_column = do_check,
			value_column = `count`,
			sample_column = sample,
			gene_column = symbol,
			pass_fit = T,
			just_discovery = T,
			approximate_posterior_inference = T,
			approximate_posterior_analysis = T,
			cores = 30,
			additional_parameters_to_save = c("intercept", "sigma_raw", "sigma_intercept", "sigma_slope", "sigma_sigma")
		)


counts_stan_MPI =

	counts %>%
	# # Subsample for test
	# inner_join(
	# 	bind_rows (
	# 		(.) %>% distinct(symbol, `Cell type category`) %>% filter(`Cell type category` == "house_keeping") %>% sample_n(200),
	# 		(.) %>% distinct(symbol) %>% sample_n(100)
	# 	) %>% distinct(symbol)
	# )%>%

	mutate(sample_idx = as.integer(sample)) %>%

	dplyr::left_join(
		(.) %>%
			distinct(symbol) %>%
			mutate( idx_MPI = head( rep(1:shards, (.) %>% nrow %>% `/` (shards) %>% ceiling ), n=(.) %>% nrow) ),
		by="symbol"
	) %>%
	arrange(idx_MPI, symbol) %>%

	# Decide start - end location
	group_by(idx_MPI) %>%
	do(
		(.) %>%
			dplyr::left_join(
				(.) %>%
					distinct(idx_MPI, sample, symbol) %>%
					arrange(idx_MPI, symbol) %>%
					count(idx_MPI, symbol) %>%
					mutate(end = cumsum(n)) %>%
					mutate(start = c(1, .$end %>% rev() %>% `[` (-1) %>% rev %>% `+` (1))),
				by=c("symbol", "idx_MPI")
			)
	) %>%
	ungroup()

data_for_stan_MPI = list(

	N = counts_stan_MPI %>% distinct(sample, symbol, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% max,
	M = counts_stan_MPI %>% distinct(start, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% max,
	G = counts_stan_MPI %>% distinct(symbol) %>% nrow(),
	S = counts_stan_MPI %>% distinct(sample) %>% nrow(),
	D = 1, #counts_stan_MPI %>% distinct(`Data base`) %>% nrow(),
	G_per_shard = counts_stan_MPI %>% distinct(symbol, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% as.array,
	n_shards = min(shards, counts_stan_MPI %>% distinct(idx_MPI) %>% nrow),
	G_per_shard_idx = c(0, counts_stan_MPI %>% distinct(symbol, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% cumsum),


	counts =
		foreach(idx = counts_stan_MPI %>% distinct(idx_MPI) %>% pull(idx_MPI), .combine = full_join) %do% {
			counts_stan_MPI %>%
				filter(idx_MPI == idx) %>%
				select(`count`) %>%
				setNames(idx) %>%
				mutate(i = 1:n())
		} %>%
		select(-i) %>%
		replace(is.na(.), 0) %>%
		mutate_if(is.numeric, as.integer) %>%
		as_matrix() %>% t,

	symbol_end =
		foreach(idx = counts_stan_MPI %>% distinct(idx_MPI) %>% pull(idx_MPI), .combine = full_join) %do% {
			counts_stan_MPI %>%
				filter(idx_MPI == idx) %>%
				distinct(end) %>%
				rbind(0, .) %>%
				#mutate(start = start - ( (.) %>% head(n=1) %>% pull(start)) + 1) %>%
				setNames(idx) %>%
				mutate(i = 1:n())
		} %>%
		select(-i) %>%
		replace(is.na(.), 0) %>%
		mutate_if(is.numeric, as.integer) %>%
		as_matrix() %>% t,

	symbols =
		foreach(idx = counts_stan_MPI %>% distinct(idx_MPI) %>% pull(idx_MPI), .combine = full_join) %do% {
			counts_stan_MPI %>%
				filter(idx_MPI == idx) %>%
				distinct(symbol) %>%
				mutate(symbol = as.character(symbol)) %>%
				setNames(idx) %>%
				mutate(i = 1:n())
		} %>%
		select(-i),

	sample_idx =
		foreach(idx = counts_stan_MPI %>% distinct(idx_MPI) %>% pull(idx_MPI), .combine = full_join) %do% {
			counts_stan_MPI %>%
				filter(idx_MPI == idx) %>%
				select(sample_idx) %>%
				setNames(idx) %>%
				mutate(i = 1:n())
		} %>%
		select(-i) %>%
		replace(is.na(.), 0) %>%
		mutate_if(is.numeric, as.integer) %>%
		as_matrix() %>% t,

	symbol_ct_idx_expanded =
		counts_stan_MPI %>%
		distinct(idx_MPI, symbol, `symbol original`, `Cell type category`) %>%
		mutate(G = 1:n()) %>% filter(`Cell type category` != "house_keeping") %>%
		mutate(C = `Cell type category` %>% as.factor %>% as.integer) %>%
		mutate(M = `symbol original` %>% as.factor %>% as.integer) %>%
		arrange(M),

	symbol_ct_idx =
		counts_stan_MPI %>%
		distinct(idx_MPI, symbol, `symbol original`, `Cell type category`) %>%
		mutate(G = 1:n()) %>% filter(`Cell type category` != "house_keeping") %>%
		mutate(C = `Cell type category` %>% as.factor %>% as.integer) %>%
		mutate(M = `symbol original` %>% as.factor %>% as.integer) %>%
		arrange(M) %>%
		select(G, C, M) %>%
		spread(C, G) %>%
		select(-M)
)

# }
# load("data_for_stan_MPI_level_1.RData")

counts_stan_MPI %>% distinct(symbol, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% as.array %>%
	paste(collapse=" ") %>%
	sprintf("Genes per shard %s", .) %>%
	print

save(list=c("data_for_stan_MPI", "counts_stan_MPI", "counts", "tree"), file=sprintf("dev/level_%s_input_%s.RData", my_level, my_run))

fileConn<-file("~/.R/Makevars")
writeLines(c("CXX14FLAGS += -O3","CXX14FLAGS += -DSTAN_THREADS", "CXX14FLAGS += -pthread"), fileConn)
close(fileConn)
Sys.setenv("STAN_NUM_THREADS" = n_cores %>% divide_by(3) %>% floor)
nb_model_MPI = stan_model("/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/RNAseq-noise-model/stan/negBinomial_tidy_MPI.stan")

Sys.time()
fit_MPI =
	sampling(
		nb_model_MPI,
		data = data_for_stan_MPI,
		chains=3, cores = 3, iter=300, warmup=200, save_warmup = FALSE,
		pars = c("lambda", "sigma_raw", "sigma_intercept", "sigma_slope", "sigma_sigma", "lambda_mu", "lambda_sigma", "exposure_rate")
	)
save(fit_MPI, file= sprintf("dev/fit_MPI_level%s_%s.RData", my_level, my_run))
Sys.time()

#load("dev/fit_MPI_level1.RData")


}

counts = readRDS(file="dev/counts_infer_NB.rds") %>% dplyr::filter(level == my_level)

counts_run = counts %>% left_join( (.) %>% distinct(symbol) %>% mutate(run = sample(1:5, n(), replace = T)))

process_fit = function(my_level, run) {

	load(sprintf("dev/fit_MPI_level%s_%s.RData", my_level,  run))
	load(sprintf("dev/level_%s_input_%s.RData", my_level,  run))

	# Parse fit and ave data
	fit_MPI %>%
		rstan::summary() %$% summary %>%
		as_tibble(rownames="par") %>%
		filter(grepl("lambda|sigma", par)) %>%
		separate(par, c("par", "G"), sep="\\[|\\]") %>%
		filter(par %in% c("lambda", "sigma_raw")) %>%
		select(par, G, "50%") %>% spread(par, "50%") %>%
		mutate(G = G %>% as.integer) %>%
		left_join(
			data_for_stan_MPI$symbols %>%
				gather(idx_MPI, symbol) %>%
				drop_na() %>% mutate(G=1:n())
		) %>%
		left_join( counts_stan_MPI %>% distinct(symbol, `Cell type category`)  ) %>%

		# Convert to gamma
		mutate(
			shape = nb_to_gamma(lambda %>% exp, 1/exp(sigma_raw)) %$% shape,
			rate  = nb_to_gamma(lambda %>% exp, 1/exp(sigma_raw)) %$% rate
		) %>%

		# Calculate confidence interval
		mutate(CI_low = qgamma(0.025, shape = shape, rate = rate)) %>%
		mutate(CI_high = qgamma(0.975, shape = shape, rate = rate)) %>%

		# QQ plots

		left_join(counts) %>%
		#do_parallel_start(n_cores, "symbol") %>%
		do({

			`%>%` = magrittr::`%>%`
			library(tidyverse)
			library(magrittr)

			(.) %>%
				group_by(symbol) %>%
				do(
					(.) %>%
						arrange(`count normalised`) %>%
						mutate(
							predicted_NB =
								qnbinom(
									ppoints(`count normalised`),
									size=.$sigma_raw %>% unique %>% exp %>% `^` (-1),
									mu=.$lambda %>% unique %>% exp
								)
						)
				) %>%
				ungroup()
		}) %>%
		#do_parallel_end() %>%
		mutate(`log of error` = (`count normalised` - predicted_NB) %>% abs %>% `+` (1) %>% log) %>%
		mutate(`error of log` = (log(`count normalised` + 1) - log(predicted_NB + 1)) ) %>%
		mutate(`error scaled` =  ((`count normalised` - predicted_NB) / (`count normalised` + 1) )) %>%
		left_join(
			(.) %>%
				group_by(symbol) %>%
				summarise(`gene error mean` = `error of log` %>% abs %>% mean)
		) %>%

		mutate(
			symbol =
				gsub(
					(.) %>%
						distinct(`Cell type category`) %>%
						pull(1) %>%
						paste("_", sep="") %>%
						paste(collapse="|"),
					"",
					symbol
				)
		) %>%

		# Normalise
		left_join(
			fit_MPI %>%
				rstan::summary() %$% summary %>%
				as_tibble(rownames="par") %>%
				filter(grepl("exposure", par)) %>%
				separate(par, c("par", "S"), sep="\\[|\\]") %>%
				select(par, S, "50%") %>%
				rename(exposure = `50%`) %>%
				select(-par) %>%
				left_join(
					counts_stan_MPI %>%
						distinct(sample, sample_idx) %>%
						rename(S = sample_idx) %>%
						mutate_if(is.integer, as.character)
				)
		) %>%

		# Recalculate `count normalised`
		mutate(`count normalised bayes` = `count` / exp(exposure)) %>%

		mutate(stan_run = !!run)

}

nb_to_gamma = function(mu, phi){
	list(
		"shape" = (mu*phi)/(mu+phi),
		"rate" = phi/(mu+phi)
	)
}

# Execute fit
calculate =
	list(
		do_inference(counts_run, 1 ),
		do_inference(counts_run, 2 ),
		do_inference(counts_run, 3 ),
		do_inference(counts_run, 4 ),
		do_inference(counts_run, 5 )
	)

options(future.globals.maxSize= (20000)*1024^2)

(1:5) %>%
	map(~ process_fit(my_level, .x)) %>%
	do.call(bind_rows, .) %>%

	mutate(level = !!my_level) %>%

	# Replace the category house_keeping
	mutate(`house keeping` = `Cell type category` == "house_keeping") %>%
	rename(temp = `Cell type category`) %>%
	left_join(
		(.) %>%
			filter(temp != "house_keeping") %>%
			distinct(sample, temp, level) %>%
			rename(`Cell type category` = temp)
	) %>%
	select(-temp) %>%
	filter(`Cell type category` %>% is.null %>% `!`) %>%

	# Keep one HKG per stan_sun
	filter((!`house keeping`) | stan_run==1) %>%

	# # Keep one copy lambda sd of house keeping genes
	# {
	# bind_rows(
	# 	(.) %>% filter(!`house keeping`),
	# 	(.) %>% filter(`house keeping`) %>%
	# 		group_by(symbol) %>%
	# 		mutate(
	# 			lambda = lambda %>% mean,
	# 			sigma_raw = sigma_raw %>% mean,
	# 			exposure = exposure %>% mean,
	# 			`gene error mean` = `gene error mean` %>% mean,
	# 		) %>%
	# 		ungroup
	# )
	# 	} %>%

		# Save
		write_csv(sprintf("dev/fit_level%s.csv", my_level))
