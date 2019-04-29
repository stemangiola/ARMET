# test new ARMET
library(tidyverse)
library(foreach)
library(magrittr)
source("R/ARMET_tc.R")

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

		# Correct house keeping formatting
		mutate(`house keeping` = `Cell type category` == "house_keeping") %>%
		rename(`Cell type category old` = `Cell type category`) %>%
		left_join(
			(.) %>%
				distinct(sample, `Cell type category old`) %>%
				filter(`Cell type category old` != "house_keeping") %>%
				rename(`Cell type category` = `Cell type category old`)
		) %>%

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

		# Re put house keeping
		mutate(`Cell type category` = ifelse(`house keeping`, "house_keeping", `Cell type category`)) %>%
		select(-`Cell type category old`) %>%

		# Eliminate n_pass info
		select(-n_pass)

}



markers =	read_csv("docs/markers.csv")

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

house_keeping =
	ref %>% filter(`Cell type category` == "house_keeping" ) %>%
	distinct(`symbol`) %>%
	pull(1) %>%
	head(n=200)


reference =

	ref %>%
	inner_join(
		(.) %>%
			distinct(sample) %>%
			filter( !grepl(sample_blacklist %>% paste(collapse="|"), sample))
	) %>%

	# Decrese number of samples
	inner_join(
		(.) %>%
			filter(level ==1) %>%
			decrease_replicates(n_pass=1) %>%
			decrease_replicates(n_pass=2) %>%
			distinct(sample)
	) %>%


	# Get house keeping and markwrs
	{
		bind_rows(
			(.) %>%	filter(	`symbol` %in% 	house_keeping 	) %>%	select(-level) %>% distinct(),
			(.) %>% inner_join(markers %>% distinct(symbol, level))
		)
	} %>%
	select(  -contains("idx")) %>%
	mutate(`read count` = `read count` %>% as.integer)

######################################
######################################

# Create mix
set.seed(123)
reps = 1
mix_source =
	ref %>%

	# Create the combinations
	group_by(level) %>%
	do(
		combn((.) %>% filter(`Cell type category` != "house_keeping") %>% distinct(`Cell type category`) %>% pull(1),m = 2)  %>%
			t %>%
			as_tibble
	) %>%
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
		library(tidyverse)
		library(magrittr)

		sample_blacklist = c("666CRI", "972UYG", "344KCP", "555QVG", "370KKZ", "511TST", "13816.11933", "13819.11936", "13817.11934", "13818.11935", "096DQV", "711SNV")
		cc = (.)
		my_ref = 	ref %>%
			distinct(sample, `Cell type category`) %>%
			filter( ! grepl(sample_blacklist %>% paste(collapse="|"), sample))

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
			left_join(ref) %>%
			distinct(`symbol`, `read count normalised bayes`, `Cell type category`, sample, level) %>%
			mutate(run = cc %>% distinct(run) %>% pull(1))
	}) %>%
	#multidplyr::collect() %>%
	ungroup() %>%

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

	# # Correct for house keeping gene catwegory
	# rename(`Cell type category OLD` = `Cell type category`) %>%
	# left_join(
	# 	(.) %>%
	# 		filter(`Cell type category OLD` != "house_keeping") %>%
	# 		distinct(sample, `Cell type category OLD`) %>%
	# 		rename(`Cell type category` = `Cell type category OLD`)
	# 	) %>%

	# Eliminate duplicated of difference levels
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
	spread(`symbol`, `read count mix`)

# Run ARMET
res = ARMET_tc(mix)
#save(res, file="temp_res.RData")

res %$%
	proportions %>%
	ungroup() %>%

	# Filter only relevant
	inner_join(	mix_source %>% distinct(level, `Cell type category`, sample_mix, combination) %>% rename(sample = sample_mix) ) %>%

	# Add real prop
	separate(combination, c("cta","ctb"), sep=" ") %>%
	rowwise %>% mutate(ct1 = min(cta, ctb), ct2 = max(cta, ctb)) %>%
	mutate(real = ifelse(`Cell type category` %in% c(ct1, ct2), 0.5, 0)) %>%
	ungroup %>%
	unite(combination, c("ct1", "ct2")) %>%

	# Plot
	{
		{ (.) %>%
			ggplot(aes(y=`.value`, x=sample, color=`Cell type category`)) +
			geom_errorbar(aes(ymin = `.lower`, ymax=`.upper`), width=0) + facet_wrap(~level + combination + real) } %>%
			ggsave(.,filename =  "level_1_test_CI.png", device = "png")
		(.)
	} %>%

	rowwise %>%
	mutate(is_outside = ifelse(real %>% between( .lower ,  .upper), 0, 1)) %>%
	mutate(error = min(abs(.lower - real), abs(.upper - real)) * is_outside) %>%
	mutate(error = ifelse(error < 1e-4, 0, error)) %>%
	{
		{(.) %>% ggplot(aes(y=error, x=combination, label=sample)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size=0.2) + facet_wrap(~ level, scales = "free_x") + my_theme} %>%
			ggsave(filename = "level_1_test_error.png", device = "png", width = 8)
		(.)
	} %>%
	group_by(combination) %>%
	summarise(`error mean` = error %>% mean)


