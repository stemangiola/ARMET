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
			distinct(`symbol`, `read count normalised bayes`, `Cell type category`, sample) %>%
			mutate(run = cc %>% distinct(run) %>% pull(1))
	}) %>%
	#multidplyr::collect() %>%
	ungroup()

mix =
	mix_source %>%
	group_by(run) %>%
	do(
		(.) %>%
			distinct(`symbol`, `read count normalised bayes`, `Cell type category`, run) %>%
			spread(`Cell type category`, `read count normalised bayes`) %>%
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
save(res, file="temp_res.RData")

res %$% proportions %>%
	ungroup() %>%

	# Add real prop
	mutate(dummy = gsub("^[0-9]+_", "", sample)) %>%
	separate(dummy, c("ct1","ct2"), remove = F, sep=" ") %>%
	mutate(ct1 = ifelse(ct1 %in% c("endothelial", "epithelial", "fibroblast"), ct1, "immune_cell")) %>%
	mutate(ct2 = ifelse(ct2 %in% c("endothelial", "epithelial", "fibroblast"), ct2, "immune_cell")) %>%
	rowwise %>% mutate(real = ifelse(`Cell type category` %in% c(ct1, ct2), 0.5, 0)) %>%
	unite(combination, c("ct1", "ct2")) %>%
	ggplot(aes(y=`.value`, x=sample, color=`Cell type category`)) + geom_errorbar(aes(ymin = `.lower`, ymax=`.upper`), width=0) + facet_wrap(~combination + real)
