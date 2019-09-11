# Create mix_base withn NA - imputed values

source(	"https://gist.githubusercontent.com/stemangiola/dd3573be22492fc03856cd2c53a755a9/raw/a7479898357de6e109419be49fe264be7775e9a9/tidy_extensions.R")
n_cores = 20

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

ToDataFrameTypeColFull = function(tree, ...){
	tree %>%
		Clone() %>%
		{
			t = (.)
			foreach(l=1:(t %$% Get("level") %>% max), .combine = bind_rows) %do% {
				data.tree::Clone(t) %>%
					{ data.tree::Prune(., function(x) x$level <= l + 1); . } %>%
					data.tree::ToDataFrameTypeCol(...) %>%
					as_tibble
			}
		} %>%
		distinct() %>%
		{ if("level_3" %in% ((.) %>% colnames)) (.) %>% mutate(level_3 = ifelse(level_3 %>% is.na, level_2, level_3)) else (.) } %>%
		{ if("level_4" %in% ((.) %>% colnames)) (.) %>% mutate(level_4 = ifelse(level_4 %>% is.na, level_3, level_4)) else (.) } %>%
		{ if("level_5" %in% ((.) %>% colnames)) (.) %>% mutate(level_5 = ifelse(level_5 %>% is.na, level_4, level_5)) else (.) } %>%
		{ if("level_6" %in% ((.) %>% colnames)) (.) %>% mutate(level_6 = ifelse(level_6 %>% is.na, level_5, level_6)) else (.) } %>%
		select(..., everything())
}


ref =
	read_csv("dev/ref_1_2_3_4.csv") %>%
	inner_join(
		(.) %>%
			distinct(sample) %>%
			filter( !grepl(sample_blacklist %>% paste(collapse="|"), sample))
	)

all_genes =
	ref %>%
	filter(level==1) %>%
	inner_join( (.) %>% distinct(symbol, `Cell type category`) %>% count(symbol) %>% filter(n == max(n)) ) %>%
	distinct(symbol) %>% pull(1)

ref2 = ref %>% filter(symbol %in% all_genes)

mix_base =
	tree %>% ToDataFrameTypeColFull() %>% select(-level_6, -level_1) %>% rename(`Cell type category` = level_5) %>%
	left_join(
		ref2 %>% distinct(`Cell type category`, level, symbol, `read count normalised bayes`, sample) %>%
			spread(symbol, `read count normalised bayes`) %>%
			gather(symbol, `read count normalised bayes`, -c(1:3))
	) %>%

	# Level 5

	# do_parallel_start(
	# 	.f = ~ {
	# 		values = .x %>% drop_na() %>% pull(`read count normalised bayes`)
	# 		if(length(values) > 0)
	#
	# 			bind_rows(
	# 				.x %>% filter(`read count normalised bayes` %>% is.na %>% `!`),
	# 				.x %>% filter(`read count normalised bayes` %>% is.na ) %>%
	# 					mutate(`read count normalised bayes` = sample(values, size = n(), replace = T))
	# 			)
	#
	# 		else (.)
	#
	# 	} ,
	# 	`Cell type category`, symbol
	# ) %>%

	do_parallel_start(n_cores, "Cell type category") %>%
	do({
		`%>%` = magrittr::`%>%`
		library(tidyverse)
		library(magrittr)

		(.) %>%
			group_by(`Cell type category`, symbol) %>%
			do({
				values = (.) %>% drop_na() %>% pull(`read count normalised bayes`)
				if(length(values) > 0)

					bind_rows(
						(.) %>% filter(`read count normalised bayes` %>% is.na %>% `!`),
						(.) %>% filter(`read count normalised bayes` %>% is.na ) %>%
							mutate(`read count normalised bayes` = sample(values, size = n(), replace = T))
					)

				else (.)

			})
	}) %>%
	do_parallel_end() %>%

	# Level 4

	# do_parallel_start(
	# 	.f = ~ {
	# 		values = .x %>% drop_na() %>% pull(`read count normalised bayes`)
	# 		if(length(values) > 0)
	#
	# 			bind_rows(
	# 				.x %>% filter(`read count normalised bayes` %>% is.na %>% `!`),
	# 				.x %>% filter(`read count normalised bayes` %>% is.na ) %>%
	# 					mutate(`read count normalised bayes` = sample(values, size = n(), replace = T))
	# 			)
	#
	# 		else (.)
	#
	# 	} ,
	# 	`level_4`, symbol
	# ) %>%


	do_parallel_start(n_cores, "level_4") %>%
	do({
		`%>%` = magrittr::`%>%`
		library(tidyverse)
		library(magrittr)

		(.) %>%
			group_by(`level_4`, symbol) %>%
			do({
				values = (.) %>% drop_na() %>% pull(`read count normalised bayes`)
				if(length(values) > 0)

					bind_rows(
						(.) %>% filter(`read count normalised bayes` %>% is.na %>% `!`),
						(.) %>% filter(`read count normalised bayes` %>% is.na ) %>%
							mutate(`read count normalised bayes` = sample(values, size = n(), replace = T))
					)

				else (.)

			})
	}) %>%
	do_parallel_end() %>%

	# Level 3

	# do_parallel_start(
	# 	.f = ~ {
	# 		values = .x %>% drop_na() %>% pull(`read count normalised bayes`)
	# 		if(length(values) > 0)
	#
	# 			bind_rows(
	# 				.x %>% filter(`read count normalised bayes` %>% is.na %>% `!`),
	# 				.x %>% filter(`read count normalised bayes` %>% is.na ) %>%
	# 					mutate(`read count normalised bayes` = sample(values, size = n(), replace = T))
	# 			)
	#
	# 		else (.)
	#
	# 	} ,
	# 	`level_3`, symbol
	# ) %>%

	do_parallel_start(n_cores, "level_3") %>%
	do({
		`%>%` = magrittr::`%>%`
		library(tidyverse)
		library(magrittr)

		(.) %>%
			group_by(`level_3`, symbol) %>%
			do({
				values = (.) %>% drop_na() %>% pull(`read count normalised bayes`)
				if(length(values) > 0)

					bind_rows(
						(.) %>% filter(`read count normalised bayes` %>% is.na %>% `!`),
						(.) %>% filter(`read count normalised bayes` %>% is.na ) %>%
							mutate(`read count normalised bayes` = sample(values, size = n(), replace = T))
					)

				else (.)

			})
	}) %>%
	do_parallel_end() %>%

	# Level 2

	# do_parallel_start(
	# 	.f = ~ {
	# 		values = .x %>% drop_na() %>% pull(`read count normalised bayes`)
	# 		if(length(values) > 0)
	#
	# 			bind_rows(
	# 				.x %>% filter(`read count normalised bayes` %>% is.na %>% `!`),
	# 				.x %>% filter(`read count normalised bayes` %>% is.na ) %>%
	# 					mutate(`read count normalised bayes` = sample(values, size = n(), replace = T))
	# 			)
	#
	# 		else (.)
	#
	# 	} ,
	# 	`level_2`, symbol
	# ) %>%

	do_parallel_start(n_cores, "level_2") %>%
	do({
		`%>%` = magrittr::`%>%`
		library(tidyverse)
		library(magrittr)

		(.) %>%
			group_by(`level_2`, symbol) %>%
			do({
				values = (.) %>% drop_na() %>% pull(`read count normalised bayes`)
				if(length(values) > 0)

					bind_rows(
						(.) %>% filter(`read count normalised bayes` %>% is.na %>% `!`),
						(.) %>% filter(`read count normalised bayes` %>% is.na ) %>%
							mutate(`read count normalised bayes` = sample(values, size = n(), replace = T))
					)

				else (.)

			})
	}) %>%
	do_parallel_end() %>%

	#attach further info
	left_join(
		ref %>%
			distinct(`symbol`, `house keeping`)
	) %>%

	# there is redundancy I don't know why
	group_by(sample, symbol, level) %>% slice(1) %>% ungroup()

saveRDS(mix_base, file="dev/mix_base.RDS")
