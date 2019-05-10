# test new ARMET
library(tidyverse)
library(foreach)
library(magrittr)


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
		axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
		strip.background = element_blank(),
		axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
		axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
	)

ref = read_csv("docs/ref.csv")

######################################
# Variable should be already set
######################################


sample_blacklist = c("666CRI", "972UYG", "344KCP", "555QVG", "370KKZ", "511TST", "13816.11933", "13819.11936", "13817.11934", "13818.11935", "096DQV", "711SNV")

decrease_replicates = function(my_df, n_pass=1){

	my_df %>%

		# Add n_pass info
		mutate(n_pass = n_pass) %>%

		# Detect redundant
		multidplyr::partition(`Cell type category`, `Data base`) %>%
		#group_by(`Cell type category`, `Data base`) %>%
		do({
			`%>%` = magrittr::`%>%`
			library(tidyverse)
			library(magrittr)
			library(foreach)
			source("https://gist.githubusercontent.com/stemangiola/dd3573be22492fc03856cd2c53a755a9/raw/e4ec6a2348efc2f62b88f10b12e70f4c6273a10a/tidy_extensions.R")
			source("https://gist.githubusercontent.com/stemangiola/90a528038b8c52b21f9cfa6bb186d583/raw/4ac3a7d74d4f7971bd8dd8eb470279eb7f42bbf8/transcription_tool_kit.R")

			(.) %>%
				get_redundant_pair(
					distance_feature = "symbol",
					redundant_feature = "sample",
					value_column = "read count normalised bayes",
					log_transform = T
				) %>%

				# For immune cells repete twice
				{
					if((.) %>% distinct(n_pass) %>% pull(1) == 1 & (.) %>% distinct(`Cell type category`) %>% pull(1) == "immune_cell" )
						(.) %>%
						get_redundant_pair(
							distance_feature = "symbol",
							redundant_feature = "sample",
							value_column = "read count normalised bayes",
							log_transform = T
						)
					else
						(.)
				}

		}) %>%

		collect %>% ungroup %>%

		# Eliminate n_pass info
		select(-n_pass)

}



#markers =	read_csv("docs/markers_pass2.csv")
markers =	read_csv("docs/markers_pass.csv")

## Trouble shoot
# plot_df =
# 	ref %>%
# 	filter(filt_for_calc) %>%
# 	decrease_replicates %>%
# 	filter(	`symbol` %in%  markers ) %>%
# 	bind_rows(
# 		mix %>%
# 			gather(symbol, `read count`, -sample) %>%
# 			mutate(`Cell type category` = "mix") %>%
# 			mutate(`read count normalised bayes` = `read count`)
# 	) %>%
# 	add_MDS_components(
# 		replicates_column = "symbol",
# 		cluster_by_column = "sample",
# 		value_column = "read count normalised bayes",
# 		log_transform = T
# 	)
#
# plot_df %>%
# 	filter(`Cell type category` != "house_keeping") %>%
# 	distinct(`PC 1`, `PC 2`, `Cell type category`, `Data base`, sample) %>%
# 	{
# 		(.) %>% ggplot(aes(
# 		x=`PC 1`, y = `PC 2`,
# 		color=`Cell type category`,
# 		label=sample, ds=`Data base`
# 		)) +
# 		geom_point() + my_theme
# 	} %>%
# 	plotly::ggplotly()


reference =

	ref %>%
	inner_join(
		(.) %>%
			distinct(sample) %>%
			filter( !grepl(sample_blacklist %>% paste(collapse="|"), sample))
	) %>%

	# Decrese number of samples
	# inner_join(
	# 	(.) %>%
	# 		filter(level ==1) %>%
	# 		decrease_replicates(n_pass=1) %>%
	# 		decrease_replicates(n_pass=2) %>%
	# 		distinct(sample)
	# ) %>%
	inner_join(
		(.) %>%
			distinct(sample, `Cell type category`) %>%
			group_by( `Cell type category`) %>%
			slice(1) %>%
			ungroup
	) %>%


	# Get house keeping and markwrs
	{
		bind_rows(
			(.) %>% inner_join( markers %>% distinct(level, symbol)),
			(.) %>% filter(`house keeping`)
		)
	} %>%

	# decrease number of house keeping
	anti_join(
		(.) %>% filter(`house keeping`) %>% distinct(symbol) %>% slice(100) #sample_frac(0.7)
	) %>%

	select(  -contains("idx")) %>%
	mutate(`read count` = `read count` %>% as.integer)

######################################
######################################
# library(parallel)
# n_cores = system("nproc", intern = TRUE) %>% as.integer
# cl = n_cores %>% makeCluster( manual=FALSE, outfile='log.txt')
# clusterEvalQ(cl, library(tidyverse))
# clusterEvalQ(cl, library(magrittr))
# clusterExport(cl, "ref")
# registerDoParallel(cl)
# multidplyr::set_default_cluster(cl)


# Create mix
set.seed(123)
reps = 1

sample_blacklist = c("666CRI", "972UYG", "344KCP", "555QVG", "370KKZ", "511TST", "13816.11933", "13819.11936", "13817.11934", "13818.11935", "096DQV", "711SNV")

my_ref = 	ref %>%
	distinct(sample, `Cell type category`) %>%
	filter( ! grepl(sample_blacklist %>% paste(collapse="|"), sample))

mix_source =
	ref %>%

	# Create the combinations
	group_by(level) %>%
	do({
		my_level = (.) %>% distinct(level) %>% pull(1)
		combn(
			(.) %>%
				distinct(`Cell type category`) %>%
				anti_join(ref %>% distinct(`Cell type category`, level) %>% filter(level < my_level ) %>% distinct(`Cell type category`)) %>%
				pull(1),
			m = 2
		)  %>%
			t %>%
			as_tibble
	}) %>%
	ungroup %>%
	mutate(`#` = 1:n()) %>%

	# Make more runs
	right_join(
		tibble(`#` = (1: ((.) %>% nrow)) %>% rep(reps))  %>%
			mutate(run = 1:n())
	) %>%
	select(-`#`) %>%

	#multidplyr::partition(run) %>%
	group_by(run) %>%
	do({
		`%>%` = magrittr::`%>%`

		cc = (.)

		bind_rows(
			my_ref %>%
				filter(`Cell type category` == (cc %>% pull(2))) %>%
				sample_n(1) %>%
				distinct(sample, `Cell type category`),
			my_ref %>%
				filter(`Cell type category` == (cc %>% pull(3))) %>%
				sample_n(1) %>%
				distinct(sample, `Cell type category`)
		) %>%
			mutate(run = cc %>% distinct(run) %>% pull(1)) %>%
			mutate(level = cc %>% distinct(level) %>% pull(1))
	}) %>%
	#multidplyr::collect() %>%
	ungroup() %>%

	# Again solving the problem with house keeping genes
	left_join(ref %>% distinct(`symbol`, `read count normalised bayes`, `Cell type category`, sample, level, `house keeping`))  %>%

	# Add mix_sample
	left_join(
		(.) %>%
		group_by(run) %>%
			do(
				(.) %>%
					distinct(`symbol`, `read count normalised bayes`, `Cell type category`, run, level) %>%
					spread(`Cell type category`, `read count normalised bayes`) %>%
					drop_na %>%
					mutate(combination = names((.))[4:5] %>% paste(collapse=" ")) %>%
					setNames(c("symbol", "run", "level", "1", "2", "combination")) %>%
					mutate( `read count mix` = ( (`1` + `2`) / 2 ) %>% as.integer ) %>%
					unite(sample_mix, c("run", "combination"), remove = F)
			) %>%
			ungroup() %>%
			distinct(symbol, run, level, combination, sample_mix,  `read count mix`)
	) %>%

	# Eliminate duplicated of difference levels for example house keepng genes
	group_by(run, symbol, sample) %>%
	arrange(level) %>%
	slice(1) %>%
	ungroup


mix =
	mix_source %>%
	distinct(symbol, sample_mix, `read count mix`) %>%
	drop_na %>%
	spread(`sample_mix`, `read count mix`) %>%
	drop_na %>%
	gather(sample_mix, `read count mix`, -symbol) %>%
	spread(`symbol`, `read count mix`) %>%
	rename(sample = sample_mix)

# Run ARMET
source("R/ARMET_tc.R")
res = ARMET_tc(mix)
#save(list = c("res", "mix_source", "mix"), file="temp_res_pass2_run2.RData")

res %$%
	proportions %>%
	ungroup() %>%

	# Filter only relevant
	#inner_join(	mix_source %>% distinct(level, `Cell type category`, sample_mix, combination) %>% rename(sample = sample_mix) ) %>%
	inner_join(	mix_source %>% distinct(level, sample_mix, combination) %>% rename(sample = sample_mix) ) %>%

	# Add real prop
	separate(combination, c("cta","ctb"), sep=" ") %>%
	rowwise %>% mutate(ct1 = min(cta, ctb), ct2 = max(cta, ctb)) %>%
	mutate(real = ifelse(`Cell type category` %in% c(ct1, ct2), 0.5, 0)) %>%
	ungroup %>%
	unite(combination, c("ct1", "ct2"),sep=" ", remove = F) %>%

	# Plot
	{
		{ (.) %>%
			ggplot(aes(y=`.value`, x=sample, color=`Cell type category`)) +
			geom_errorbar(aes(ymin = `.lower`, ymax=`.upper`), width=0) + facet_wrap(~level + combination + real) +
				theme(strip.text.x = element_text(size = 5)) } %>%
			ggsave(.,filename =  "level_1_test_CI.png", device = "png")
		(.)
	} %>%

	rowwise %>%
	mutate(is_outside = ifelse(real %>% between( .lower ,  .upper), 0, 1)) %>%
	mutate(error = min(abs(.lower - real), abs(.upper - real)) * is_outside) %>%

	# Analysis error on worst pairs

	group_by(sample, real) %>%
	arrange(error %>% desc) %>%
	slice(1) %>%
	ungroup %>%
	select(sample, `Cell type category`, real, error) %>%
	group_by(sample) %>%
	mutate(`error mean` = error %>% mean) %>%
	select(-error) %>%
	spread(real, `Cell type category`) %>%
	arrange(`error mean` %>% desc) %>%
	unite(pair, c("0", "0.5"), sep=" ") %>%
	ungroup %>%
	{
			((.) %>% ggplot(aes(x=pair, y=`error mean`)) + geom_boxplot() + geom_jitter() + my_theme) %>%
				ggsave(filename = "level_1_test_error.png", device = "png", width = 8)
		(.)
	} %>%

	# mutate(error = ifelse(error < 1e-4, 0, error)) %>%
	# {
	# 	{(.) %>% ggplot(aes(y=error, x=combination, label=sample)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size=0.2) + facet_wrap(~ level + real, scales = "free_x") + my_theme} %>%
	# 		ggsave(filename = "level_1_test_error.png", device = "png", width = 8)
	# 	(.)
	# } %>%


	# Correct temporary mistake
	# anti_join(ref %>% filter(level <2) %>% distinct(`Cell type category`) %>% setNames("ct1") %>% mutate(level=2)) %>%
	# anti_join(ref %>% filter(level <2) %>% distinct(`Cell type category`) %>% setNames("ct2") %>% mutate(level=2)) %>%
	separate(pair, c("ct1", "ct2"), remove = F, sep=" ") %>%
	mutate(`error min` = 0.05) %>% # `error mean` %>% min) %>%
	mutate(`error mean relative` = (`error mean` / `error min`) %>% sqrt) %>%
	left_join(
		read_csv("docs/num_markers_based_on_error_levels_1_2_first_run.csv") %>%
			distinct(ct1, ct2, `n markers`) %>% rename(`n markers pass 1` = `n markers`)
	) %>%
	mutate(`n markers` = (`error mean relative` * `n markers pass 1` ) %>% ceiling) %>%
	mutate(`n markers` = max(20, `n markers`)) %>%
	write_csv("docs/num_markers_based_on_error_levels_1_2_second_run.csv")



