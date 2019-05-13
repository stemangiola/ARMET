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

# Setup table of name conversion
level_df = foreach( l = list(
	c("b_cell", "immune_cell", "b_cell"),
	c("b_memory", "immune_cell", "b_cell", "b_memory"),
	c("b_naive", "immune_cell", "b_cell", "b_naive"),
	c("dendritic",  "immune_cell", "dendritic"),
	c("dendritic_m",  "immune_cell", "dendritic"),
	c("dendritic_m_immature",  "immune_cell", "dendritic", "dendritic_m_immature"),
	c("dendritic_m_mature",  "immune_cell", "dendritic", "dendritic_m_mature"),
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
	c("t_CD8", "immune_cell","t_cell"),
	c("t_CD8_memory_effector", "immune_cell","t_cell", "t_CD8_memory_effector"),
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

set.seed(123)

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
					mutate(pair = names((.))[4:5] %>% paste(collapse=" ")) %>%
					setNames(c("symbol", "run", "level", "1", "2", "pair")) %>%
					mutate( `read count mix` = ( (`1` + `2`) / 2 ) %>% as.integer ) %>%
					unite(sample_mix, c("run", "pair"), remove = F)
			) %>%
			ungroup() %>%
			distinct(symbol, run, level, pair, sample_mix,  `read count mix`)
	) %>%

	# Eliminate duplicated of difference levels for example house keepng genes
	group_by(run, symbol, sample) %>%
	arrange(level) %>%
	slice(1) %>%
	ungroup

get_input_data = function(markers, reps, pass){

	list(
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
					#slice(1) %>%
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
				(.) %>% filter(`house keeping`) %>% distinct(symbol) %>% slice(500) #sample_frac(0.7)
			) %>%

			select(  -contains("idx")) %>%
			mutate(`read count` = `read count` %>% as.integer)
		,
		mix =
			mix_source %>%
			distinct(symbol, sample_mix, `read count mix`) %>%
			drop_na %>%
			spread(`sample_mix`, `read count mix`) %>%
			drop_na %>%
			gather(sample_mix, `read count mix`, -symbol) %>%
			spread(`symbol`, `read count mix`) %>%
			rename(sample = sample_mix)
	) %>%
		{
		 input =	(.)
		 save(input, file=sprintf("docs/input_test_pass_%s.RData", pass))
		 (.)
		}


}

load("fit_df_feature_selection.RData")

marker_df =
	get_marker_df_source(fit_df_1, 1, 0.5) %>%
	rbind(get_marker_df_source(fit_df_2, 2, 0.4)) %>%
	separate(pair, c("ct1", "ct2"), sep=" ", remove = F)

get_markers_number = function(pass, res, num_markers_previous_level){

	if(pass == 0 | (num_markers_previous_level %>% is.null))

		marker_df %>% distinct(pair, ct1, ct2, level) %>% mutate(`n markers` = 20)

	else

	res %$%
		proportions %>%
		ungroup() %>%

		# Filter only relevant, leave 3rd party comparisons
		inner_join(	mix_source %>% distinct(level, sample_mix, pair) %>% rename(sample = sample_mix) ) %>%

		# Add real prop
		separate(pair, c("cta","ctb"), sep=" ") %>%
		rowwise %>% mutate(ct1 = min(cta, ctb), ct2 = max(cta, ctb)) %>%
		mutate(real = ifelse(`Cell type category` %in% c(ct1, ct2), 0.5, 0)) %>%
		ungroup %>%
		unite(pair, c("ct1", "ct2"),sep=" ", remove = F) %>%

		# Plot
		{
			{ (.) %>%
					ggplot(aes(y=`.value`, x=sample, color=`Cell type category`)) +
					geom_errorbar(aes(ymin = `.lower`, ymax=`.upper`), width=0) + facet_wrap(~level + pair + real) +
					theme(strip.text.x = element_text(size = 5)) } %>%
				ggsave(.,filename = sprintf("pass_%s_test_CI.png", pass), device = "png")
			(.)
		} %>%

		rowwise %>%
		mutate(is_outside = ifelse(real %>% between( .lower ,  .upper), 0, 1)) %>%
		mutate(error = min(abs(.lower - real), abs(.upper - real)) * is_outside) %>%
		ungroup() %>%

		# Penalise divergent
		mutate(error = ifelse(!converged, 1, error)) %>%

		{
			tbl = (.)

			tbl1 = tbl %>%
				group_by(sample, level) %>%
				arrange(error %>% desc) %>%
				slice(1:2) %>%
				select(sample, `Cell type category`, real, `error`, level) %>%
				summarise(`error mean` = error %>% mean, pair = paste(`Cell type category`, collapse=" ")) %>%
				ungroup() %>%
				separate(pair, c("ct1", "ct2"), sep=" ", remove = F) %>%
				bind_rows(
					(.) %>% mutate(dummy = ct1, ct1=ct2, ct2=dummy) %>% select(-dummy) %>% unite(pair, c("ct1", "ct2"), sep=" ", remove = F)
				) %>%
				inner_join(marker_df %>% distinct(level, ct1, ct2))

				bind_rows(
					tbl1,
					tbl %>%
						bind_rows(
							(.) %>% mutate(dummy = ct1, ct1=ct2, ct2=dummy) %>% select(-dummy) %>% unite(pair, c("ct1", "ct2"), sep=" ", remove = F)
						) %>%
						anti_join( tbl1 %>% distinct(pair)) %>%
						group_by(pair) %>%
						arrange(error %>% desc) %>%
						slice(1) %>%
						rename(`error mean` = error) %>%
						ungroup
				)
		} %>%

		group_by(pair, level) %>%
		summarise(`error mean` = `error mean` %>% mean) %>%
		{
			((.) %>% ggplot(aes(x=pair, y=`error mean`)) + geom_boxplot() + geom_jitter() + my_theme) %>%
				ggsave(filename = sprintf("pass_%s_test_error.png", pass), device = "png", width = 8)
			(.)
		} %>%
		ungroup() %>%

		mutate(`error min` = 0.05) %>% # `error mean` %>% min) %>%
		mutate(`error mean relative` = (`error mean` / `error min`) %>% sqrt) %>%
		left_join(
			num_markers_previous_level %>%
				distinct(pair, ct1, ct2, `n markers`) %>% rename(`n markers pass 1` = `n markers`)
		) %>%
		mutate(`n markers` = (`error mean relative` * `n markers pass 1` ) %>% ceiling) %>%
		mutate(`n markers` = ifelse(`n markers` < 20, 20, `n markers`)) %>%
		{
			(.) %>% write_csv(sprintf("docs/num_markers_based_on_error_levels_1_2_pass_%s.csv", pass))
			(.)
		}
}


get_marker_df_source = function(fit_df, level, fit_threshold, lambda_threshold =4){

	fit_df %>%
		left_join(level_df %>% mutate(`level 0` = "root")) %>%
		filter(`Cell type category` != "house_keeping") %>%
		inner_join(
			(.) %>%
				select(!!sprintf("level %s", level -1), !!sprintf("level %s", level)) %>%
				distinct %>%
				group_by_at(1) %>%
				summarise(n = n()) %>%
				filter(n > 1)
		) %>%
		group_by(!!sym(sprintf("level %s", level -1))) %>%
		do(
			gtools::permutations(
				n=(.) %>% select(!!sprintf("level %s", level)) %>% distinct %>% drop_na %>% nrow,
				r=2,
				v=(.) %>% select(!!sprintf("level %s", level)) %>% distinct %>% drop_na %>% pull(1),
				repeats.allowed=F
			) %>%
				as_tibble
		) %>%
		mutate(comparison = 1:n()) %>%
		unite(pair, c("V1", "V2"), remove = F, sep=" ") %>%
		gather(which, `Cell type category`, -!!sprintf("level %s", level -1), -comparison, -pair) %>%
		left_join(
			fit_df %>% distinct(symbol, `Cell type category`, CI_low, CI_high)
		) %>%
		mutate(relevant_CI = ifelse(which == "V1", CI_low, CI_high)) %>%
		group_by(comparison, pair, symbol) %>%
		arrange(which) %>%
		summarise(delta = diff(relevant_CI)) %>%
		separate(pair, c("Cell type category", "other Cell type category"), sep=" ", remove = F) %>%
		left_join( fit_df %>% distinct(`Cell type category`, symbol, lambda, `gene error mean`) ) %>%
		filter(lambda > !!lambda_threshold) %>%
		filter(`gene error mean` < fit_threshold) %>%
		arrange(delta) %>%
		mutate(rank = 1:n()) %>%
		ungroup() %>%
		mutate(level = !!level)
}

get_markers_df = function(markers_number, pass){
	marker_df %>%

		separate(pair, c("ct1", "ct2"), remove = F, sep=" ") %>%
		left_join(markers_number, by=c("ct1", "ct2", "level") ) %>%

		# Filter markers
		filter(rank < `n markers`) %>%

		# Filter markers in upper levels that have been selected for lower levels
		inner_join(
			(.) %>%
				distinct(symbol, level) %>%
				group_by(symbol) %>%
				arrange(level %>% desc) %>%
				slice(1) %>%
				ungroup
		) %>%

		# Write table
		{
			(.)  %>%	write_csv(sprintf("docs/markers_pass%s.csv", pass))
			(.)
		}
}

source("R/ARMET_tc.R")
reps = 1
##################################
# Pass 0
##################################

n_markers_0 = get_markers_number(0, NULL, NULL)
res_0 =
	n_markers_0 %>%
	get_markers_df(0) %>%
	get_input_data(reps = reps, pass = 0) %>%
	{	ARMET_tc((.)$mix, (.)$reference) }

##################################
# Pass 1
##################################

n_markers_1 = get_markers_number(1, res_0, n_markers_0)
res_1 =
	n_markers_1 %>%
	get_markers_df(1) %>%
	get_input_data(reps = reps, pass = 1) %>%
	{	ARMET_tc((.)$mix, (.)$reference) }

##################################
# Pass 2
##################################

n_markers_2 = get_markers_number(2, res_1, n_markers_1)
res_2 =
	n_markers_2 %>%
	get_markers_df(2) %>%
	get_input_data(reps = reps, pass = 2) %>%
	{	ARMET_tc((.)$mix, (.)$reference) }

