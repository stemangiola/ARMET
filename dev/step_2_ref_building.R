library(tidyverse)
library(magrittr)
library(purrr)
library(furrr)
library(data.tree)
library(foreach)
library(ARMET)
library(SIBERG)
source("~/PostDoc/ppcSeq/R/do_parallel.R")

ifelse_pipe = function(.x, .p, .f1, .f2 = NULL) {
	switch(.p %>% `!` %>% sum(1),
				 as_mapper(.f1)(.x),
				 if (.f2 %>% is.null %>% `!`)
				 	as_mapper(.f2)(.x)
				 else
				 	.x)

}

options(future.globals.maxSize = 50000 * 1024 ^ 2)

give_rank_to_ref = function(fit_df, level, tree, lambda_threshold = 5) {

	n_ct = fit_df %>% distinct(`Cell type category`) %>% nrow

	fit_df %>%
		left_join(
			tree %>%
			Clone() %>% ARMET::ToDataFrameTypeColFull(fill = T, "name") %>% as_tibble() %>% rename(`Cell type formatted` = name)
		) %>%
		filter(`Cell type category` != "house_keeping") %>%

		# Select only branches that have childs
		inner_join(
			(.) %>%
				select(
					!!sprintf("level_%s", level),
					!!sprintf("level_%s", level + 1)
				) %>%
				distinct %>%
				group_by_at(1) %>%
				summarise(n = n()) %>%
				filter(n > 1),
			by=sprintf("level_%s", level)
		) %>%

		# For each category
		group_by(!!sym(sprintf("level_%s", level))) %>%
		do(
			gtools::permutations(
				n = (.) %>% select(!!sprintf("level_%s", level + 1)) %>% distinct %>% drop_na %>% nrow,
				r = 2,
				v = (.) %>% select(!!sprintf("level_%s", level + 1)) %>% distinct %>% drop_na %>% pull(1),
				repeats.allowed = F
			) %>%
				as_tibble
		) %>%
		mutate(comparison = 1:n()) %>%
		unite(pair, c("V1", "V2"), remove = F, sep = " ") %>%
		gather(which,
					 `Cell type category`,
					 -!!sprintf("level_%s", level),
					 -comparison,
					 -pair) %>%
		left_join(fit_df %>% distinct(symbol, `Cell type category`, CI_low, CI_high)) %>%

		# Temporary for cell nameserror

		filter(symbol %>% is.na %>% `!`) %>%
		anti_join(
			(.) %>% distinct(pair, `Cell type category`, symbol) %>% count(pair, symbol) %>% arrange(n) %>% filter(n ==
																																																						 	1) %>% ungroup %>% distinct(pair)
		) %>%

		mutate(relevant_CI = ifelse(which == "V1", CI_low, CI_high)) %>%
		group_by(comparison, pair, symbol) %>%
		arrange(which) %>%
		summarise(delta = diff(relevant_CI)) %>%
		separate(
			pair,
			c("Cell type category", "other Cell type category"),
			sep = " ",
			remove = F
		) %>%
		left_join(
			fit_df %>% distinct(
				`Cell type category`,
				symbol,
				lambda_log,
				sigma_inv_log,
				`gene error mean`,
				regression,
				bimodality_NB, bimodality_NB_diff, `marker too noisy`, `symbol too noisy`
			)
		) %>%
		filter(lambda_log > !!lambda_threshold) %>%
		# {
		# 	if( (.) %>% filter(`Cell type category` == "eosinophil") %>% nrow %>% `>` (0)) browser()
		# 	(.)
		# } %>%

		# If a marker exists for more cell types (can happen) none of them can be noisy otherwise we screw up the whole gene
		ungroup() %>%

		# For markers Filter if the marker itself is noisy
		group_by(symbol) %>%
		add_tally() %>%
		mutate(`high trans noisy frac` = sum(`marker too noisy`)/n) %>%
		ungroup() %>%
		filter(`high trans noisy frac` < 0.1666667) %>%

		# For non markers Filter symbols that have more than 10% of noisy markers for any cell type
		filter(`symbol too noisy` <= (n_ct*0.20) %>% floor ) %>%

		# Filter markers that have competitor genes way to highly transcribed (e.g., eosinophil 1000 neutrophil 10, but monocytes 100000)
		# group_by(symbol) %>%
		# mutate(`mean lambda_log` = lambda_log %>% mean) %>%
		# ungroup() %>%
		# filter(lambda_log > `mean lambda_log`) %>%


		group_by(comparison, pair) %>%
		arrange(delta) %>%
		mutate(rank = 1:n()) %>%
		ungroup() %>%
		mutate(level = !!level)
}

create_ref = function(ref, markers) {
	ref %>%

		# Get house keeping and markwrs
		left_join(markers %>% distinct(level, symbol, ct1, ct2, rank) %>% filter(rank < 500)) %>%
		filter(`house keeping` | rank %>% is.na %>% `!`) %>%

		# Filter out symbol if present both in markers and house keeping (strange but happens)
		anti_join((.) %>% filter(`house keeping` & (ct1 %>% is.na %>% `!`)) %>% distinct(symbol)) %>%

		select(-contains("idx")) %>%
		mutate(`count` = `count` %>% as.integer)
}

plot_trends = function(input.df, symbols) {
	input.df %>%
		filter(symbol %in% symbols) %>%
		filter(level == 1) %>%
		#filter(`Cell type category` == "endothelial") %>%
		ggplot(
			aes(
				x = `predicted_NB` + 1,
				y = `count normalised bayes` + 1,
				label = symbol,
				color = `regression`
			)
		) +
		geom_abline(intercept = 0,
								slope = 1,
								color = "grey") +
		geom_jitter() +
		geom_smooth(method = "lm", formula = y ~ x + I(x ^ 2)) +
		facet_wrap( ~ symbol + `Cell type category`, scales = "free")  +
		expand_limits(y = 1, x = 1) +
		scale_y_log10() +
		scale_x_log10()

}

mixture_mu_ratios = function(x){
	ss = SIBERG::SIBER(x, model='NB')
	max(ss[1], ss[2]) / (min(ss[1], ss[2]) + 1)
}

regression_coefficient = function(x){

	mdf = x %>%
		get_NB_qq_values %>%
		mutate(a = `count normalised bayes` %>% `+` (1) %>% log,
					 b = predicted_NB %>% `+` (1) %>% log)

	lm(a ~  b + I(b ^ 2), data = mdf)$coefficients[3]
}


ref =
	list(
		read_csv("dev/fit_level1.csv"),
		read_csv("dev/fit_level2.csv"),
		read_csv("dev/fit_level3.csv")
	) %>%
	map_dfr( ~ .x) %>%

	# Bug after I deleted FANTOM5 I have to rerun infer NB. Some genes are not in all cell types anynore
	# Other bug
	filter(`Cell type category` %>% is.na %>% `!`) %>%
	group_by(level) %>%
	do(
		(.) %>% inner_join( (.) %>% distinct(symbol, `Cell type category`) %>% count(symbol) %>% filter(n == max(n)), by="symbol" )
	) %>%
	ungroup() %>%

	# Eliminate house keeping that are not in all levels
	inner_join(
		bind_rows(
			(.) %>% filter(!`house keeping`),
			(.) %>% filter(`house keeping`) %>%
				distinct(symbol, level) %>%
				count(symbol) %>%
				filter(n == max(n))
		) %>%
			distinct(symbol),
		by="symbol"
	)

ref_2 =
	ref %>%

	# Add regression coefficients
	do_parallel_start(20, "symbol") %>%
	do({

		`%>%` = magrittr::`%>%`
		library(tidyverse)
		library(magrittr)

		get_NB_qq_values = function(input.df) {
			input.df = input.df %>% arrange(`count normalised bayes`)

			my_predicted_NB =
				qnbinom(
					# If 1 sample, just use median
					switch(
						input.df %>% nrow %>% `>` (1) %>% `!` %>% sum(1),
						ppoints(input.df$`count normalised bayes`),
						0.5
					),
					size = input.df$sigma_inv_log %>% unique %>% exp %>% `^` (-1),
					mu = input.df$lambda_log %>% unique %>% exp
				)

			input.df %>%	mutate(predicted_NB = my_predicted_NB)
		}

		(.) %>%
			group_by(`Cell type category`, symbol, level) %>%
			do({
				mdf = (.) %>%
					get_NB_qq_values %>%
					mutate(a = `count normalised bayes` %>% `+` (1) %>% log,
								 b = predicted_NB %>% `+` (1) %>% log)

				mdf %>%
					mutate(regression = lm(a ~  b + I(b ^ 2), data = (.))$coefficients[3])

			})
	}) %>%
	do_parallel_end() %>%

	# Calculate bimodality
	# mutate(bimodality = modes::bimodality_coefficient(`count normalised bayes` )) %>%
	# ungroup %>%

	# Calculate bimodality
	do_parallel_start(20, "symbol") %>%
	do({

		`%>%` = magrittr::`%>%`
		library(tidyverse)
		library(magrittr)

		siberg_iterative = function(x) {

			if(x %>% unique %>% length %>% `<` (5)) return(c(NA, NA))



			mu = NA
			max_i = ceiling(length(x)/10)
			i = 0
			while (mu %>% is.na | i <= max_i) {
				res = SIBERG::SIBER(x, model='NB')

				BI = res[7]
				mu = res[1]
				x = x[-1]
				i = i+1

			}


			if(mu %>% is.na & x %>% length %>% `<` (8)) return(c(NA, NA))

			return(c(
				max(res[1], res[2]) / (min(res[1], res[2]) + 1),
				res[7]
			))
		}



		(.) %>%
			group_by(`Cell type category`, symbol, level) %>%
			do({
				res = siberg_iterative((.) %>% pull(`count normalised bayes`) %>% as.integer)
				(.) %>% mutate(bimodality_NB = res[2], bimodality_NB_diff = res[1])

			})
	}) %>%
	do_parallel_end() %>%

	# Classify noisy genes
	mutate(
		`marker too noisy` =

			# Regression too parabolic
			(regression >= 0.2 & sigma_inv_log >= 0.5) |

			# Too bimodal
			( (bimodality_NB > 0.8 & bimodality_NB_diff > 20) | bimodality_NB_diff > 100)

	) %>%

	# If NA is bad
	mutate(`marker too noisy` = ifelse(`marker too noisy` %>% is.na, T, `marker too noisy`)) %>%

	# Salve all really lowly transcribed
	mutate(`marker too noisy` = ifelse(lambda_log <= log(10), F, `marker too noisy`)) %>%
	mutate(`marker too noisy` = ifelse(sigma_inv_log <= 0, F, `marker too noisy`)) %>%

	# Summarise by each gene
	left_join(
		(.) %>%
			distinct(level, `Cell type category`, symbol, `marker too noisy`) %>%
			group_by(level, symbol) %>%
			summarise(`symbol too noisy` = sum(`marker too noisy`))
	)

saveRDS(ref_2, file='dev/ref_2.rds')

# Some plots
# ref_2 %>%
# 	distinct(
# 		level,
# 		symbol,
# 		`Cell type category`,
# 		regression,
# 		`gene error mean`,
# 		lambda_log,
# 		sigma_inv_log,
# 		sample
# 	) %>%
# 	count(level,
# 				symbol,
# 				`Cell type category`,
# 				regression,
# 				`gene error mean`,
# 				lambda_log,
# 				sigma_inv_log) %>%
# 	unite(symbol_ct, c("symbol", "Cell type category"), remove = F) %>%
# 	arrange(`regression` %>% desc) %>%
# 	mutate(symbol_ct = factor(symbol_ct, unique(symbol_ct))) %>%
#
# 	mutate(quarantine = regression >= 0.2 & sigma_inv_log >= 0.5) %>%
#
# 	sample_frac(0.01) %>%
# 	filter(level == 1) %>%
# 	{
# 		ggplot((.),
# 					 aes(
# 					 	x = `sigma_inv_log`,
# 					 	y = `regression` ,
# 					 	color = `quarantine`,
# 					 	size = n,
# 					 	label = symbol
# 					 )) +
# 			geom_point(alpha = 0.4, width = 25) +
# 			facet_wrap( ~ `Cell type category`) +
# 			theme(
# 				axis.title.x = element_blank(),
# 				axis.text.x = element_blank(),
# 				axis.ticks.x = element_blank()
# 			)
# 	} %>%
# 	plotly::ggplotly()
#
#
# ref_2 %>% plot_trends("MAP6")
# ref_2 %>% plot_boxplot("MAP6")

my_tree = ARMET::create_tree_object(ref_2)

ref_3 =
	create_ref(
		ref_2,
		give_rank_to_ref(ref_2 %>% filter(level == 1), 1, my_tree) %>%
			rbind(give_rank_to_ref(ref_2 %>% filter(level == 2), 2, my_tree)) %>%
			rbind(give_rank_to_ref(ref_2 %>% filter(level == 3), 3, my_tree)) %>%
			#rbind(give_rank_to_ref(ref_2 %>% filter(level == 4), 4, my_tree)) %>%
			separate(pair, c("ct1", "ct2"), sep = " ", remove = F)
	)

# (ref_3 %>% filter(ct1 == "epithelial" & ct2 == "endothelial" & level ==1) %>%
# 		arrange(rank) %>% inner_join( (.) %>% distinct(symbol) %>% slice(1:50)) %>%
# 		ggplot(aes(x=`Cell type category`, y=`count normalised bayes`+1, color=`marker too noisy`)) +
# 		geom_jitter() + facet_wrap(~symbol ) + scale_y_log10() ) %>% plotly::ggplotly()

ARMET_ref = ref_3 %>% mutate_if(is.character, as.factor) %>% mutate(`count normalised bayes` = `count normalised bayes` %>% as.integer)

# Select correct columns
columns = c(  "level"     ,                  "sample"    ,                  "symbol"    ,
	  "Cell type formatted"    ,     "count"   ,               "count normalised bayes",
	  "Data base"          ,         "lambda_log"        ,              "sigma_inv_log"    ,
	          "CI_low"   ,
	  "CI_high"     ,                "gene error mean"     ,        "house keeping"    ,
	  "Cell type category"  ,               "predicted_NB" ,
	  "a"            ,               "b"                  ,         "regression"  ,
	  "bimodality_NB"   ,            "bimodality_NB_diff"  ,        "marker too noisy"    ,
	  "symbol too noisy"   ,         "ct1"             ,            "ct2"  ,
	  "rank" )
ARMET_ref = ARMET_ref %>% select(!!columns)

# # Adjust too little noise to avoid overfitting
# ARMET_ref =
# 	ARMET_ref %>%
# 	mutate(rigma_raw_intercept = 1.52, sigma_inv_log_slope = -0.4, sigma_inv_log_sigma = 1.2) %>%
# 	mutate(sigma_inv_log_regressed = ( rigma_raw_intercept + lambda_log * sigma_inv_log_slope) ) %>%
# 	mutate(sigma_inv_log_minimum =  sigma_inv_log_regressed - sigma_inv_log_sigma)

save(ARMET_ref, file="data/ARMET_ref.RData", compress = "xz")


# ARMET_ref %>%
# 	distinct(level, symbol, `Cell type category`, lambda_log, sigma_inv_log) %>%
# 	ggplot(aes(x = lambda_log, y = sigma_inv_log, color=`Cell type category`)) +
# 	geom_point() + facet_wrap(~level) +
# 	geom_smooth(method="lm") +
# 	geom_abline(intercept = 1.52, slope=-0.4) +
# 	geom_abline(intercept = 1.52 - 1.2, slope=-0.4)
#
# N52_ARMET_T$signatures %>%
# 	map_dfr(~ .x) %>%
# 	filter(!query & !`house keeping`) %>%
# 	filter(level ==1) %>%
# 	distinct(symbol, `Cell type category`, lambda_log, sigma_inv_log) %>%
# 	stan_glm(sigma_inv_log ~ lambda_log, data = .)
#
# N52_ARMET_T$signatures %>%
# 	map_dfr(~ .x) %>%
# 	filter(!query & !`house keeping`) %>%
# 	filter(level ==1) %>%
# 	distinct(symbol, `Cell type category`, lambda_log, sigma_inv_log) %>%
# 	brm(data = .,
# 			bf(family = student,
# 		sigma_inv_log ~ lambda_log,
# 		nu = 4
# 			)
# 	)
