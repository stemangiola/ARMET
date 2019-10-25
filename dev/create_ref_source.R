library(tidyverse)
library(magrittr)
library(purrr)
library(furrr)
library(data.tree)
library(foreach)
library(ARMET)
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

give_rank_to_ref = function(fit_df, level, lambda_threshold = 5) {
	ToDataFrameTypeColFull = function(tree, ...) {
		tree %>%
			Clone() %>%
			{
				t = (.)
				foreach(l = 1:(t %$% Get("level") %>% max), .combine = bind_rows) %do% {
					data.tree::Clone(t) %>%
						{
							data.tree::Prune(., function(x)
								x$level <= l + 1)
							.
						} %>%
						data.tree::ToDataFrameTypeCol(...) %>%
						as_tibble
				}
			} %>%
			distinct() %>%
			{
				if ("level_3" %in% ((.) %>% colnames))
					(.) %>% mutate(level_3 = ifelse(level_3 %>% is.na, level_2, level_3))
				else
					(.)
			} %>%
			{
				if ("level_4" %in% ((.) %>% colnames))
					(.) %>% mutate(level_4 = ifelse(level_4 %>% is.na, level_3, level_4))
				else
					(.)
			} %>%
			{
				if ("level_5" %in% ((.) %>% colnames))
					(.) %>% mutate(level_5 = ifelse(level_5 %>% is.na, level_4, level_5))
				else
					(.)
			} %>%
			{
				if ("level_6" %in% ((.) %>% colnames))
					(.) %>% mutate(level_6 = ifelse(level_6 %>% is.na, level_5, level_6))
				else
					(.)
			} %>%
			select(..., everything())
	}

	n_ct = fit_df %>% distinct(`Cell type category`) %>% nrow

	fit_df %>%
		left_join(
			Clone(ARMET::tree) %>% ToDataFrameTypeColFull("name") %>% as_tibble() %>% rename(`Cell type formatted` = name)
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
				filter(n > 1)
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
				lambda,
				sigma_raw,
				`gene error mean`,
				regression,
				bimodality_NB, bimodality_NB_diff, `marker too noisy`, `symbol too noisy`
			)
		) %>%
		filter(lambda > !!lambda_threshold) %>%
		# {
		# 	if( (.) %>% filter(`Cell type category` == "eosinophil") %>% nrow %>% `>` (0)) browser()
		# 	(.)
		# } %>%

		# If a marker exists for more cell types (can happen) none of them can be noisy otherwise we screw up the whole gene
		ungroup() %>%

		# Filter symbols that have more than 10% of noisy markers for any cell type
		filter(`symbol too noisy` <= (n_ct*0.15) %>% floor ) %>%

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
		anti_join((.) %>% filter(`house keeping` &
														 	(ct1 %>% is.na %>% `!`)) %>% distinct(symbol)) %>%

		select(-contains("idx")) %>%
		mutate(`read count` = `read count` %>% as.integer)
}

plot_trends = function(input.df, symbols) {
	input.df %>%
		filter(symbol %in% symbols) %>%
		filter(level == 1) %>%
		#filter(`Cell type category` == "endothelial") %>%
		ggplot(
			aes(
				x = `predicted_NB` + 1,
				y = `read count normalised bayes` + 1,
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

plot_boxplot = function(input.df, symbols) {
	input.df %>%
		filter(symbol %in% symbols) %>%
		filter(level == 1) %>%
		#filter(`Cell type category` == "endothelial")  %>%
		ggplot(
			aes(
				x = `Cell type category`,
				y = `read count normalised bayes` + 1,
				label = symbol,
				color = `regression`
			)
		) +
		geom_jitter() +
		geom_boxplot() +
		facet_wrap( ~ symbol + `Cell type category`, scales = "free")  +
		expand_limits(y = 1, x = 1) +
		scale_y_log10()
}

mixture_mu_ratios = function(x){
	ss = SIBERG::SIBER(x, model='NB')
	max(ss[1], ss[2]) / (min(ss[1], ss[2]) + 1)
}

regression_coefficient = function(x){

	mdf = x %>%
		get_NB_qq_values %>%
		mutate(a = `read count normalised bayes` %>% `+` (1) %>% log,
					 b = predicted_NB %>% `+` (1) %>% log)

	lm(a ~  b + I(b ^ 2), data = mdf)$coefficients[3]
}

sample_blacklist = c(
	"666CRI",
	"972UYG",
	"344KCP",
	"555QVG",
	"370KKZ",
	"511TST",
	"13816.11933",
	"13819.11936",
	"13817.11934",
	"13818.11935",
	"096DQV",
	"711SNV",
	"counts.Ciliary%20Epithelial%20Cells%2c%20donor3.CNhs12009.11399-118D4"   ,
	"counts.Iris%20Pigment%20Epithelial%20Cells%2c%20donor1.CNhs12596.11530-119I9",
	"ENCFF890DJO"
)

ref =
	read_csv("dev/ref_1_2_3_4.csv") %>%
	inner_join((.) %>%
						 	distinct(sample) %>%
						 	filter(!grepl(
						 		sample_blacklist %>% paste(collapse = "|"), sample
						 	))) %>%
	filter(`Data base` != "FANTOM5") %>%

	# Filter dendritic that look too much like monocytes
	filter(!(
		(`Cell type formatted` == "dendritic_myeloid" & `Data base` == "Immune Singapoor") |
			(`Cell type formatted` == "dendritic_myeloid" & `Data base` == "bloodRNA")
	))

ref_2 =
	ref %>%

	# Add regression coefficients
	do_parallel_start(20, "symbol") %>%
	do({

		`%>%` = magrittr::`%>%`
		library(tidyverse)
		library(magrittr)

		get_NB_qq_values = function(input.df) {
			input.df = input.df %>% arrange(`read count normalised bayes`)

			predicted_NB =
				qnbinom(
					# If 1 sample, just use median
					switch(
						input.df %>% nrow %>% `>` (1) %>% `!` %>% sum(1),
						ppoints(input.df$`read count normalised bayes`),
						0.5
					),
					size = input.df$sigma_raw %>% unique %>% exp %>% `^` (-1),
					mu = input.df$lambda %>% unique %>% exp
				)

			input.df %>%	mutate(predicted_NB = predicted_NB)
		}

		(.) %>%
			group_by(`Cell type category`, symbol, level) %>%
			do({
				mdf = (.) %>%
					get_NB_qq_values %>%
					mutate(a = `read count normalised bayes` %>% `+` (1) %>% log,
								 b = predicted_NB %>% `+` (1) %>% log)

				mdf %>%
					mutate(regression = lm(a ~  b + I(b ^ 2), data = (.))$coefficients[3])

			})
	}) %>%
	do_parallel_end() %>%

	# Calculate bimodality
	# mutate(bimodality = modes::bimodality_coefficient(`read count normalised bayes` )) %>%
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
				res = siberg_iterative((.) %>% pull(`read count normalised bayes`) %>% as.integer)
				(.) %>% mutate(bimodality_NB = res[2], bimodality_NB_diff = res[1])

			})
	}) %>%
	do_parallel_end() %>%

	# Classify noisy genes
	mutate(
		`marker too noisy` =

			# Regression too parabolic
			(regression >= 0.2 & sigma_raw >= 0.5) |

			# Too bimodal
			( (bimodality_NB > 0.8 & bimodality_NB_diff > 20) | bimodality_NB_diff > 100)

	) %>%

	# If NA is bad
	mutate(`marker too noisy` = ifelse(`marker too noisy` %>% is.na, T, `marker too noisy`)) %>%

	# Salve all really lowly transcribed
	mutate(`marker too noisy` = ifelse(lambda <= log(10), F, `marker too noisy`)) %>%
	mutate(`marker too noisy` = ifelse(sigma_raw <= 0, F, `marker too noisy`)) %>%

	# Summarise by each gene
	left_join(
		(.) %>%
			distinct(level, `Cell type category`, symbol, `marker too noisy`) %>%
			group_by(level, symbol) %>%
			summarise(`symbol too noisy` = sum(`marker too noisy`))
	)

saveRDS(ref_2, file='ref_2.rds')

# Some plots
# ref_2 %>%
# 	distinct(
# 		level,
# 		symbol,
# 		`Cell type category`,
# 		regression,
# 		`gene error mean`,
# 		lambda,
# 		sigma_raw,
# 		sample
# 	) %>%
# 	count(level,
# 				symbol,
# 				`Cell type category`,
# 				regression,
# 				`gene error mean`,
# 				lambda,
# 				sigma_raw) %>%
# 	unite(symbol_ct, c("symbol", "Cell type category"), remove = F) %>%
# 	arrange(`regression` %>% desc) %>%
# 	mutate(symbol_ct = factor(symbol_ct, unique(symbol_ct))) %>%
#
# 	mutate(quarantine = regression >= 0.2 & sigma_raw >= 0.5) %>%
#
# 	sample_frac(0.01) %>%
# 	filter(level == 1) %>%
# 	{
# 		ggplot((.),
# 					 aes(
# 					 	x = `sigma_raw`,
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

ref_3 =
	create_ref(
		ref_2,
		give_rank_to_ref(ref_2 %>% filter(level == 1), 1) %>%
			rbind(give_rank_to_ref(ref_2 %>% filter(level == 2), 2)) %>%
			rbind(give_rank_to_ref(ref_2 %>% filter(level == 3), 3)) %>%
			rbind(give_rank_to_ref(ref_2 %>% filter(level == 4), 4)) %>%
			separate(pair, c("ct1", "ct2"), sep = " ", remove = F)
	)

# (ref_3 %>% filter(ct1 == "epithelial" & ct2 == "endothelial" & level ==1) %>%
# 		arrange(rank) %>% inner_join( (.) %>% distinct(symbol) %>% slice(1:50)) %>%
# 		ggplot(aes(x=`Cell type category`, y=`read count normalised bayes`+1, color=`marker too noisy`)) +
# 		geom_jitter() + facet_wrap(~symbol ) + scale_y_log10() ) %>% plotly::ggplotly()

ARMET_ref = ref_3 %>% mutate_if(is.character, as.factor) %>% mutate(`read count normalised bayes` = `read count normalised bayes` %>% as.integer)

save(ARMET_ref, file="data/ARMET_ref.RData", compress = "xz")
