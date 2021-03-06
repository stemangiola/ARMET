# Create mix_base withn NA - imputed values

library(tidyverse)
library(magrittr)
library(purrr)
library(furrr)
library(data.tree)
library(foreach)
library(ARMET)
#source("~/PhD/deconvolution/ARMET/R/utils.R")
#source("~/PostDoc/ppcSeq/R/do_parallel.R")
n_cores = 20

do_parallel_start = function(df, cores, partition_by){

  # Only if cores > 1
  if(cores > 1)		cl <- multidplyr::new_cluster(cores)

  df %>%
    dplyr::left_join(
      (.) %>%
        dplyr::select(!!partition_by) %>%
        dplyr::distinct() %>%
        dplyr::mutate(
          `.part` = 1:n() %>%
            magrittr::divide_by(length((.))) %>%
            magrittr::multiply_by(!!cores) %>%
            ceiling
        )
    )  %>%
    group_by(.part) %>%

    # Only if cores > 1
    ifelse_pipe(cores > 1,	~ .x %>% multidplyr::partition(cl))

}


do_parallel_end = function(.){
  (.) %>%
    # Only if cores > 1
    ifelse_pipe((.) %>% class %>% magrittr::equals("multidplyr_party_df") %>% any,	~ .x %>% dplyr::collect()) %>%
    dplyr::ungroup() %>%

    # Only if cores > 1
    dplyr::select(-`.part`)
}

ifelse_pipe = function(.x, .p, .f1, .f2 = NULL) {
	switch(.p %>% `!` %>% sum(1),
				 as_mapper(.f1)(.x),
				 if (.f2 %>% is.null %>% `!`)
				 	as_mapper(.f2)(.x)
				 else
				 	.x)

}

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

ref =	ARMET::ARMET_ref

all_genes =
	ref %>%
	filter(level==1) %>%
	inner_join(
		(.) %>% 
			distinct(symbol, `Cell type category`) %>% 
			dplyr::count(symbol) %>%
			filter(n == max(n)) 
	) %>%
	distinct(symbol) %>% pull(1)

mix_base_sparse = ref %>% filter(symbol %in% all_genes)

impute_missing = function(ref2, tree){

	yaml::yaml.load_file("data/tree.yaml") %>%
		data.tree::as.Node() %>% ToDataFrameTypeColFull() %>% select(-level_6, -level_1) %>% rename(`Cell type category` = level_5) %>%
		left_join(
			ref2 %>% distinct(`Cell type category`, level, symbol, `count scaled bayes`, sample) %>%
				spread(symbol, `count scaled bayes`) %>%
				gather(symbol, `count scaled bayes`, -c(1:3))
		) %>%

	do_parallel_start(n_cores, "Cell type category") %>%
		do({
			`%>%` = magrittr::`%>%`
			library(tidyverse)
			library(magrittr)

			(.) %>%
				group_by(`Cell type category`, symbol) %>%
				do({
					values = (.) %>% drop_na() %>% pull(`count scaled bayes`)
					if(length(values) > 0)

						bind_rows(
							(.) %>% filter(`count scaled bayes` %>% is.na %>% `!`),
							(.) %>% filter(`count scaled bayes` %>% is.na ) %>%
								mutate(`count scaled bayes` = sample(values, size = n(), replace = T))
						)

					else (.)

				})
		}) %>%
		do_parallel_end() %>%



	do_parallel_start(n_cores, "level_4") %>%
		do({
			`%>%` = magrittr::`%>%`
			library(tidyverse)
			library(magrittr)

			(.) %>%
				group_by(`level_4`, symbol) %>%
				do({
					values = (.) %>% drop_na() %>% pull(`count scaled bayes`)
					if(length(values) > 0)

						bind_rows(
							(.) %>% filter(`count scaled bayes` %>% is.na %>% `!`),
							(.) %>% filter(`count scaled bayes` %>% is.na ) %>%
								mutate(`count scaled bayes` = sample(values, size = n(), replace = T))
						)

					else (.)

				})
		}) %>%
		do_parallel_end() %>%

		# Level 3

		# do_parallel_start(
		# 	.f = ~ {
		# 		values = .x %>% drop_na() %>% pull(`count scaled bayes`)
		# 		if(length(values) > 0)
		#
		# 			bind_rows(
		# 				.x %>% filter(`count scaled bayes` %>% is.na %>% `!`),
		# 				.x %>% filter(`count scaled bayes` %>% is.na ) %>%
	# 					mutate(`count scaled bayes` = sample(values, size = n(), replace = T))
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
					values = (.) %>% drop_na() %>% pull(`count scaled bayes`)
					if(length(values) > 0)

						bind_rows(
							(.) %>% filter(`count scaled bayes` %>% is.na %>% `!`),
							(.) %>% filter(`count scaled bayes` %>% is.na ) %>%
								mutate(`count scaled bayes` = sample(values, size = n(), replace = T))
						)

					else (.)

				})
		}) %>%
		do_parallel_end() %>%

		# Level 2

		# do_parallel_start(
		# 	.f = ~ {
		# 		values = .x %>% drop_na() %>% pull(`count scaled bayes`)
		# 		if(length(values) > 0)
		#
		# 			bind_rows(
		# 				.x %>% filter(`count scaled bayes` %>% is.na %>% `!`),
		# 				.x %>% filter(`count scaled bayes` %>% is.na ) %>%
	# 					mutate(`count scaled bayes` = sample(values, size = n(), replace = T))
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
					values = (.) %>% drop_na() %>% pull(`count scaled bayes`)
					if(length(values) > 0)

						bind_rows(
							(.) %>% filter(`count scaled bayes` %>% is.na %>% `!`),
							(.) %>% filter(`count scaled bayes` %>% is.na ) %>%
								mutate(`count scaled bayes` = sample(values, size = n(), replace = T))
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
		group_by(sample, symbol, level) %>% dplyr::slice(1) %>% ungroup()

}


mix_base = impute_missing(mix_base_sparse, tree)
saveRDS(mix_base, file="dev/mix_base.RDS")


# Create noiless mix

mix_base_noiseless = impute_missing(
	mix_base_sparse %>%
		distinct(
			`Cell type category`,
			level,
			symbol,
			`count scaled bayes` = exp( lambda_log),
			sample = `Cell type category`
		),
	tree
)
saveRDS(mix_base_noiseless, file="dev/mix_base_noiseless.RDS", compress = "gzip")
