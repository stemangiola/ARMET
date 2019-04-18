# test new ARMET
library(tidyverse)
library(foreach)
load("docs/ref.RData")

# Create mix
set.seed(123)
mix_source =
	combn(ref %>% filter(`Cell type category` != "house_keeping") %>% distinct(`Cell type category`) %>% pull(1),m = 2) %>%
	t %>%
	as_tibble %>%
	mutate(`#` = 1:n()) %>%

	# Make more runs
	right_join(
		tibble(`#` = (1: ((.) %>% nrow)) %>% rep(50))  %>%
			mutate(run = 1:n())
	) %>%
	select(-`#`) %>%

	#multidplyr::partition(run) %>%
	group_by(run) %>%
	do({
		`%>%` = magrittr::`%>%`
		library(tidyverse)
		library(magrittr)

		sample_blacklist = c("666CRI", "972UYG", "344KCP", "555QVG", "370KKZ", "511TST", "13816.11933", "13819.11936", "13817.11934", "13818.11935", "096DQV", "711SNV")
		cc = (.)
		my_ref = 	ref %>%
			distinct(sample, `Cell type category`) %>%
			filter( ! grepl(sample_blacklist %>% paste(collapse="|"), sample))

		bind_rows(
			my_ref %>%
				filter(`Cell type category` == (cc %>% pull(1))) %>%
				sample_n(1) %>%
				distinct(sample),
			my_ref %>%
				filter(`Cell type category` == (cc %>% pull(2))) %>%
				sample_n(1) %>%
				distinct(sample)
		) %>%
			left_join(ref) %>%
			distinct(`symbol`, `read count normalised bayes`, `Cell type formatted`, sample) %>%
			mutate(run = cc %>% distinct(run) %>% pull(1))
	}) %>%
	#multidplyr::collect() %>%
	ungroup()

mix =
	mix_source %>%
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

# Run ARMET
source("R/ARMET_tc.R")
res = ARMET_tc(mix)
