# Create mix_base withn NA - imputed values

ref = read_csv("dev/ref_1_2_3_4.csv")

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
