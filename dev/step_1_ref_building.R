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
library(multidplyr)
library(data.tree)
library(ttBulk)
library(ARMET)
library(furrr)
library(purrr)

# Get level
args <- commandArgs(TRUE)
my_level = args[1] %>% as.integer

# # Registering parallel framework
# system("nproc", intern = TRUE) %>%
# 	as.integer %>%
# 	sprintf("Working on %s cores", .) %>%
# 	print
#
# n_cores = system("nproc", intern = TRUE) %>% as.integer %>% sum(-4)

# # Tools
# source("../../../PostDoc/ppcSeq/R/do_parallel.R")
# source("https://gist.githubusercontent.com/stemangiola/9d2ba5d599b7ac80404c753cdee04a01/raw/26e5b48fde0cd4f5b0fd7cbf2fde6081a5f63e7f/tidy_data_tree.R")

# MPI
do_inference = function(counts, my_level, my_run, approximate_posterior = F){

	counts =
		counts %>%
		filter(`Cell type category` == "house_keeping" | run == my_run) %>%
		mutate(symbol = `symbol original`)

	counts_fit =
		counts %>%

		# inner_join(
		# 	bind_rows(
		# 		(.) %>% distinct(symbol) %>% sample_n(100),
		# 		(.) %>% filter(`Cell type category` == "house_keeping") %>% distinct(symbol) %>% slice(1:100)
		# 	) %>%
		# 		distinct
		#
		# ) %>%

		select(sample, symbol, count, `Cell type category`, level) %>%
		rename(`read count` = count) %>%

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

		ARMET::ref_intercept_only(my_level, cores = ceiling(24/5), approximate_posterior = approximate_posterior)

	counts_fit %>%
		left_join(counts) %>%
		saveRDS(sprintf("dev/lambda_sigma_MPI_level_%s_%s.RData", my_level,  my_run))

	counts_fit %>% attr("fit") %>%
		saveRDS(		sprintf("dev/fit_MPI_level_%s_%s.RData", my_level,  my_run)		)

}

process_fit = function(my_level, run) {


	# Parse fit and ave data
	readRDS(sprintf("dev/lambda_sigma_MPI_level_%s_%s.RData", my_level,  run)) %>%

		# Convert to gamma
		mutate(
			shape = nb_to_gamma(lambda_log %>% exp, 1/exp(sigma_inv_log)) %$% shape,
			rate  = nb_to_gamma(lambda_log %>% exp, 1/exp(sigma_inv_log)) %$% rate
		) %>%

		# Calculate confidence interval
		mutate(CI_low = qgamma(0.025, shape = shape, rate = rate)) %>%
		mutate(CI_high = qgamma(0.975, shape = shape, rate = rate)) %>%

		# QQ plots
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
									size=.$sigma_inv_log %>% unique %>% exp %>% `^` (-1),
									mu=.$lambda_log %>% unique %>% exp
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

		# Normalise
		mutate(`count normalised bayes` = `count` / exp(exposure)) %>%

		mutate(stan_run = !!run)

}

nb_to_gamma = function(mu, phi){
	list(
		"shape" = (mu*phi)/(mu+phi),
		"rate" = phi/(mu+phi)
	)
}

counts_run =
	readRDS(file="dev/counts_infer_NB.rds") %>%
	left_join( (.) %>% distinct(`symbol original`) %>% mutate(run = sample(1:5, n(), replace = T)))


# Execute fit
plan(multiprocess)
options(future.globals.maxSize= (20000)*1024^2)

# Fit
calculate =
	1:5 %>%
	future_map(
		~ do_inference(	counts_run,	my_level,	.x ,	T	)
	)

# Process fit
(1:5) %>%
	future_map_dfr(~ process_fit(my_level, .x)) %>%

	# Add level information
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

	# Save
	write_csv(sprintf("dev/fit_level%s.csv", my_level))
