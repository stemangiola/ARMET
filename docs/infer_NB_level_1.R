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
library(parallel)
library(doParallel)
library(multidplyr)

# library(future)
# plan(multiprocess)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Get level
args <- commandArgs(TRUE)
my_level = args[1]

# Registering parallel framework
system("nproc", intern = TRUE) %>%
	as.integer %>%
	sprintf("Working on %s cores", .) %>%
	print

n_cores = system("nproc", intern = TRUE) %>% as.integer

cl = n_cores %>% makeCluster( manual=FALSE, outfile='log.txt')

clusterEvalQ(cl, library(tidyverse))
clusterEvalQ(cl, library(magrittr))

registerDoParallel(cl)
set_default_cluster(cl)

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
source("https://gist.githubusercontent.com/stemangiola/90a528038b8c52b21f9cfa6bb186d583/raw/dbd92c49fb03fb05ab0b465704b99c0a39e654d5/transcription_tool_kit.R")
source("https://gist.githubusercontent.com/stemangiola/dd3573be22492fc03856cd2c53a755a9/raw/e4ec6a2348efc2f62b88f10b12e70f4c6273a10a/tidy_extensions.R")

# Setup table of name conversion
level_df = foreach( l = list(
	c("b_cell", "immune_cell", "b_cell"),
	c("b_memory", "immune_cell", "b_cell", "b_memory"),
	c("b_naive", "immune_cell", "b_cell", "b_naive"),
	c("dendritic",  "immune_cell", "dendritic_myeloid"),
	c("dendritic_plasmacytoid",  "immune_cell", "dendritic_plasmacytoid"),
	c("dendritic_myeloid",  "immune_cell", "dendritic_myeloid"),
	c("dendritic_myeloid_immature",  "immune_cell", "dendritic_myeloid", "dendritic_myeloid_immature"),
	c("dendritic_myeloid_mature",  "immune_cell", "dendritic_myeloid", "dendritic_myeloid_mature"),
	c("endothelial", "endothelial", "endothelial", "endothelial"),
	c("eosinophil", "immune_cell", "granulocyte", "eosinophil"),
	c("epithelial", "epithelial", "epithelial", "epithelial"),
	c("fibroblast","fibroblast", "fibroblast", "fibroblast"),
	c("immune_cell", "immune_cell"),
	c("macrophage", "immune_cell", "mono_derived"),
	c("macrophage_M1", "immune_cell", "mono_derived", "macrophage_M1"),
	c("macrophage_M2", "immune_cell", "mono_derived", "macrophage_M2"),
	c("mast_cell", "immune_cell", "mast_cell", "mast_cell"),
	c("monocyte","immune_cell", "mono_derived",  "monocyte"),
	c("myeloid", "immune_cell"),
	c("natural_killer", "immune_cell", "natural_killer", "natural_killer"),
	c("neutrophil", "immune_cell", "granulocyte", "neutrophil"),
	c("t_CD4", "immune_cell","t_cell"),
	c("t_CD8", "immune_cell","t_cell", "t_CD8"),
	c("t_CD8_memory_effector", "immune_cell","t_cell", "t_CD8"),
	c("t_cell", "immune_cell","t_cell"),
	c("t_gamma_delta", "immune_cell","t_cell", "t_gamma_delta"),
	c("t_helper", "immune_cell","t_cell", "t_helper"),
	c("t_memory_central", "immune_cell","t_cell", "t_memory_central"),
	c("t_memory", "immune_cell","t_cell"),
	c("t_memory_effector", "immune_cell","t_cell", "t_memory_effector"),
	c("t_naive", "immune_cell","t_cell", "t_naive"),
	c("t_reg", "immune_cell","t_cell", "t_reg"),
	c("t_reg_memory", "immune_cell","t_cell", "t_reg")
), .combine = bind_rows) %do% {
	l %>%
		as.data.frame %>% t %>% as_tibble(.name_repair = "minimal") %>%
		setNames(c("Cell type formatted", sprintf("level %s",  1:((.)%>%ncol()-1) ) ))
} %>%
	mutate(`level 0` = "root")

ct_to_correlation_threshold =
	level_df %>%
	select(-`Cell type formatted`) %>%
	gather(label, `Cell type category`) %>%
	separate(label, c("dummy", "level")) %>%
	mutate(level = level %>% as.integer) %>%
	filter(level %in% 1:3) %>%
	group_by(`Cell type category`) %>%
	arrange(level) %>%
	slice(1) %>%
	left_join(
		tibble(level=1:6, threshold=c(0.9, 0.95, 0.99, 0.99, 0.99, 0.99))
	)

# if(0){
# Get data
counts =
	foreach(
		f = dir(
			path = "/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/RNAseq-noise-model/big_data/tibble_cellType_files",
			pattern="cellTypes_",
			full.names = T
		),
		.combine = bind_rows
	) %dopar% {
		read_csv(f) %>%
			mutate_at(vars(one_of('sample')), as.character) %>%
			mutate_at(vars(one_of('isoform')), as.character)
	} %>%

	filter(symbol %>% is.na %>% `!`) %>%
	filter(`Cell type formatted` %>% is.na %>% `!`) %>%

	# Setup Cell type category names
	left_join(
		level_df %>%
			select(`Cell type formatted`, contains(sprintf("level %s", my_level))) %>%
			setNames(c((.) %>% colnames  %>% `[` (1), "Cell type category"))

	) %>%
	filter(`Cell type category` %>% is.na %>% `!`) %>%
	mutate(level = !!my_level) %>%

	# Eliminate genes that are not in all cell types
	inner_join( (.) %>% distinct(symbol, `Cell type category`) %>% count(symbol) %>% filter(n == max(n)) ) %>%

	# Setup house keeping genes
	mutate(`Cell type category` = ifelse(symbol %in% (read_csv("docs/hk_600.txt", col_names = FALSE) %>% pull(1)), "house_keeping", `Cell type category`)) %>%

	# Median redundant
	do_parallel_start(n_cores, "symbol") %>%
	do({
		`%>%` = magrittr::`%>%`
		library(tidyverse)
		library(magrittr)

		(.) %>%
			group_by(sample, symbol, `Cell type`, `Cell type category`, `Cell type formatted`,  `Data base`) %>%
			summarise(`read count` = `read count` %>% median(na.rm = T)) %>%
			ungroup()
	}) %>%
	do_parallel_end() %>%

	# Normalise
	norm_RNAseq(
		sample_column = "sample",
		gene_column = "symbol",
		value_column = "read count",
		cpm_threshold = 0.5
	) %>%
	mutate(`read count normalised` = `read count normalised` %>% as.integer) %>%
	mutate(`read count normalised log` = `read count normalised` %>% `+` (1) %>% log) %>%

	# mutate symbol
	mutate(`symbol original` = symbol) %>%
	unite(symbol, c("Cell type category", "symbol"), remove = F) %>%

	# Mark the bimodal distributions
	do_parallel_start(n_cores, "symbol") %>%
	do({
		`%>%` = magrittr::`%>%`
		library(tidyverse)
		library(magrittr)

		(.) %>%
			group_by(symbol) %>%
			do(
				(.) %>%
					mutate(
						`bimodal p-value` =
							(.) %>%
							pull(`read count normalised log`) %>%
							diptest::dip.test() %$%
							`p.value`,

						`bimodal coefficient` =
							(.) %>%
							pull(`read count normalised log`) %>%
							modes::bimodality_coefficient(),

						`anova p-value` =
							ifelse(
								(.) %>% distinct(`Data base`) %>% nrow > 1 & (.) %>% distinct(`Cell type formatted`) %>% nrow > 1 ,
								(.) %>% aov(`read count normalised log` ~ `Cell type formatted` + `Data base`, .) %>% anova %$% `Pr(>F)` %>% `[` (1),
								ifelse(
									(.) %>% distinct(`Data base`) %>% nrow > 1,
									(.) %>% aov(`read count normalised log` ~ `Data base`, .) %>% anova %$% `Pr(>F)` %>% `[` (1),
									ifelse(
										(.) %>% distinct(`Cell type formatted`) %>% nrow > 1,
										(.) %>% aov(`read count normalised log` ~ `Cell type formatted`, .) %>% anova %$% `Pr(>F)` %>% `[` (1),
										NA
									)
								)

							)
					)
			) %>%
			ungroup()

	}) %>%
	do_parallel_end() %>%
	mutate(
		`anova p-value` = ifelse(`anova p-value` == "NaN", 1, `anova p-value`),
		`bimodal coefficient` = ifelse(`bimodal coefficient` == "NaN", 0, `bimodal coefficient`)
	) %>%
	mutate(
		`hard bimodality` =
			(`bimodal p-value` < 0.05) +
			(`bimodal coefficient` > 0.6666667) +
			(`anova p-value` < 0.05 ) >= 2
	) %>%
	mutate(`soft bimodality` = `anova p-value` < 0.0001) %>%

	# Remove redundant samples
	{

		nr =
			(.) %>%
			filter(`Cell type category` != "house_keeping") %>%
			group_by(`Cell type category`) %>%
			do({

				threshold = (.) %>% distinct(`Cell type category`) %>% left_join( ct_to_correlation_threshold ) %>% pull(threshold)


				# Remove redundant samples
				(.) %>%
					anti_join(
					(.) %>%
						distinct(symbol, sample, `read count`) %>%
						spread(sample, `read count`) %>%
						drop_na %>%
						gather(sample, `read count`, -symbol) %>%
						rename(rc = `read count`) %>%
						mutate_if(is.factor, as.character) %>%
						widyr::pairwise_cor(sample, symbol, rc, sort=T, diag = FALSE, upper = F) %>%
						filter(correlation > threshold) %>%
						distinct(item1) %>%
						rename(sample = item1)
				) %>%

				# Sample homogeneous population
				mutate( `threshold contribution` = (.) %>%  distinct(sample, `Cell type formatted`) %>% count(`Cell type formatted`) %>% pull(n) %>% quantile(0.8) %>% floor) %>%
				group_by(`Cell type formatted`) %>%
				inner_join( (.) %>% distinct(sample, `threshold contribution`) %>%  filter(row_number(sample) <= `threshold contribution`)) %>%
				ungroup() %>%
				select(- `threshold contribution`)

			}) %>%
			ungroup

		(.) %>%
			filter(`Cell type category` == "house_keeping") %>%
			inner_join( nr %>% distinct(sample)) %>%
			bind_rows( nr )

	} %>%

	mutate(symbol = symbol %>% as.factor) %>%
	mutate(sample = sample %>% as.factor)

# MPI

shards = n_cores %>% divide_by(3) %>% floor %>% multiply_by(4)
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
				select(`read count`) %>%
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

save(list=c("data_for_stan_MPI", "counts_stan_MPI", "counts", "level_df"), file=sprintf("docs/level_%s_input.RData", my_level))

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
		chains=3, iter=300, warmup=200, save_warmup = FALSE,
		pars = c("lambda", "sigma_raw", "sigma_intercept", "sigma_slope", "sigma_sigma", "lambda_mu", "lambda_sigma", "exposure_rate")
	)
save(fit_MPI, file= sprintf("docs/fit_MPI_level%s.RData", my_level))
Sys.time()

#load("docs/fit_MPI_level1.RData")

# Parse fit and ave data
fit_MPI %>%
summary() %$% summary %>%
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

# QQ plots

left_join(counts) %>%
do_parallel_start(n_cores, "symbol") %>%
do({

	`%>%` = magrittr::`%>%`
	library(tidyverse)
	library(magrittr)

	(.) %>%
		group_by(symbol) %>%
		do(
			(.) %>%
				arrange(`read count normalised`) %>%
				mutate(
					predicted_NB =
						qnbinom(
							ppoints(`read count normalised`),
							size=.$sigma_raw %>% unique %>% exp %>% `^` (-1),
							mu=.$lambda %>% unique %>% exp
						)
				)
		) %>%
		ungroup()
}) %>%
do_parallel_end() %>%
mutate(`log of error` = (`read count normalised` - predicted_NB) %>% abs %>% `+` (1) %>% log) %>%
mutate(`error of log` = (log(`read count normalised` + 1) - log(predicted_NB + 1)) ) %>%
mutate(`error scaled` =  ((`read count normalised` - predicted_NB) / (`read count normalised` + 1) )) %>%
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
	write_csv(sprintf("docs/fit_level%s.csv", my_level))
