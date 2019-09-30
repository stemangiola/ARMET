# test new ARMET
library(tidyverse)
library(foreach)
library(magrittr)
library(foreach)
library(doParallel)
registerDoParallel()
library(ARMET)
library(data.tree)

#options(error = quote({dump.frames(to.file=TRUE); q()}), show.error.locations = TRUE)

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


# out_dir = "temp"

# Get level
args <- commandArgs(TRUE)

n_markers = args[1] %>% as.integer
reps = 10
is_full_bayesian = args[2] %>% as.integer %>% as.logical
is_full_bayesian = T

levels = 1:(args[3] %>% as.integer)

#out_dir = Sys.time() %>% format("%a_%b_%d_%X") %>% gsub("[: ]", "_", .) %>% sprintf("dev/feature_selection_%s", .)
out_dir = args[2]  %>% sprintf("dev/feature_selection_%s", .)
out_dir %>% dir.create()

iterations = 500

# Variable should be already set

get_markers_number = function(pass, res, num_markers_previous_level, min_n_samples=20){

	if(pass == 0 | (num_markers_previous_level %>% is.null))

		marker_df %>% distinct(pair, ct1, ct2, level) %>% left_join( tibble(level=c(1,2, 3, 4), `n markers` = c( (min_n_samples * 2.5) %>% ceiling, min_n_samples, min_n_samples, min_n_samples)) )

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
				ggsave(.,filename = sprintf("%s/pass_%s_test_CI.png", out_dir, pass), device = "png")
			(.)
		} %>%

		rowwise %>%
		mutate(is_outside = ifelse(real %>% between( .lower ,  .upper), 0, 1)) %>%
		mutate(error = min(abs(.lower - real), abs(.upper - real)) * is_outside) %>%
		ungroup() %>%

		# Penalise divergent
		#mutate(error = ifelse(!converged, 1, error)) %>%

		# Get pairwise errors
		{

			tbl = (.)

			tbl1 = tbl %>%
				group_by(sample, level) %>%
				do(
					bind_rows(
						(.) %>%
							arrange(error %>% desc) %>%
							slice(1:2),
						(.) %>%
							filter(`real` == 0.5)
					)
				) %>%
				distinct() %>%
				select(sample, `Cell type category`, real, `error`, level) %>%
				# Take pairwise error means
				do({
					tbl = (.)

					gtools::permutations(
						n=tbl %>% nrow,
						r=2,
						v=tbl %>% pull(`Cell type category`) %>% as.character(),
						repeats.allowed=F
					) %>%
						as_tibble %>%
						rowwise %>%
						do({
							pair = (.) %>% as.character
							tbl %>%
								filter(`Cell type category` %in% pair) %>%
								summarise(`error mean` = error %>% mean, pair = paste(`Cell type category`, collapse=" "))
						})
				}) %>%
				ungroup() %>%
				separate(pair, c("ct1", "ct2"), sep=" ", remove = F) %>%

				# Add flipper combinations
				# bind_rows(
				# 	(.) %>% mutate(dummy = ct1, ct1=ct2, ct2=dummy) %>% select(-dummy) %>% unite(pair, c("ct1", "ct2"), sep=" ", remove = F)
				# ) %>%
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

		# group_by(pair, level) %>%
		# summarise(`error mean` = `error mean` %>% mean) %>%
		{
			#browser()
			((.) %>% ggplot(aes(x=pair, y=`error mean`)) + geom_boxplot() + geom_jitter() + my_theme) %>%
				ggsave(filename = sprintf("%s/pass_%s_test_error.png", out_dir, pass), device = "png", width = 8)
			(.)
		} %>%
		ungroup() %>%

		mutate(`error min` = 0.05) %>% # `error mean` %>% min) %>%
		mutate(`error mean relative` = (`error mean` / `error min`) %>% sqrt) %>%
		left_join(
			num_markers_previous_level %>%
				distinct(pair, ct1, ct2, `n markers`) %>% rename(`n markers pass 1` = `n markers`)
		) %>%
		group_by(pair, ct1, ct2, level, `n markers pass 1`) %>%
		summarise(`error mean relative mean` = `error mean relative` %>% mean) %>%
		ungroup() %>%

		mutate(`n markers` = (`error mean relative mean` * `n markers pass 1` ) %>% ceiling) %>%
		mutate(`n markers` = ifelse(`n markers` < mean(`n markers pass 1`, !!min_n_samples), mean(`n markers pass 1`, !!min_n_samples) , `n markers`)) %>%
		{
			(.) %>% write_csv(sprintf("%s/num_markers_based_on_error_levels_1_2_pass_%s.csv", out_dir, pass))
			(.)
		}
}

give_rank_to_ref = function(fit_df, level, fit_threshold, lambda_threshold =5){

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


	fit_df %>%
		left_join(  Clone(ARMET::tree) %>% ToDataFrameTypeColFull("name") %>% as_tibble() %>% rename(`Cell type formatted`= name)) %>%
		filter(`Cell type category` != "house_keeping") %>%

		# Select only branches that have childs
		inner_join(
			(.) %>%
				select(!!sprintf("level_%s", level), !!sprintf("level_%s", level + 1)) %>%
				distinct %>%
				group_by_at(1) %>%
				summarise(n = n()) %>%
				filter(n > 1)
		) %>%

		# For each category
		group_by(!!sym(sprintf("level_%s", level))) %>%
		do(
			gtools::permutations(
				n=(.) %>% select(!!sprintf("level_%s", level + 1)) %>% distinct %>% drop_na %>% nrow,
				r=2,
				v=(.) %>% select(!!sprintf("level_%s", level + 1)) %>% distinct %>% drop_na %>% pull(1),
				repeats.allowed=F
			) %>%
				as_tibble
		) %>%
		mutate(comparison = 1:n()) %>%
		unite(pair, c("V1", "V2"), remove = F, sep=" ") %>%
		gather(which, `Cell type category`, -!!sprintf("level_%s", level), -comparison, -pair) %>%
		left_join(
			fit_df %>% distinct(symbol, `Cell type category`, CI_low, CI_high)
		) %>%

		# Temporary for cell nameserror

		filter(symbol %>% is.na %>% `!`) %>%
		anti_join(
			(.) %>% distinct(pair, `Cell type category`, symbol) %>% count(pair, symbol) %>% arrange(n) %>% filter(n==1) %>% ungroup %>% distinct(pair)
		) %>%

		mutate(relevant_CI = ifelse(which == "V1", CI_low, CI_high)) %>%
		group_by(comparison, pair, symbol) %>%
		arrange(which) %>%
		summarise(delta = diff(relevant_CI)) %>%
		separate(pair, c("Cell type category", "other Cell type category"), sep=" ", remove = F) %>%
		left_join( fit_df %>% distinct(`Cell type category`, symbol, lambda, `gene error mean`) ) %>%
		filter(lambda > !!lambda_threshold) %>%
		# {
		# 	if( (.) %>% filter(`Cell type category` == "eosinophil") %>% nrow %>% `>` (0)) browser()
		# 	(.)
		# } %>%

		# If a marker exists for more cell types (can happen) none of them can be noisy otherwise we screw up the whole gene
		ungroup() %>%
		group_by(symbol) %>%
		mutate(`too noisy` = `gene error mean` %>% max %>% `>` (fit_threshold)) %>%
		ungroup %>%
		filter(!`too noisy`) %>%
		group_by(comparison, pair) %>%

		#filter(`gene error mean` < fit_threshold) %>%
		arrange(delta) %>%
		mutate(rank = 1:n()) %>%
		ungroup() %>%
		mutate(level = !!level)
}

get_markers_df = function(markers_number, pass){
	marker_df %>%

		separate(pair, c("ct1", "ct2"), remove = F, sep=" ") %>%
		left_join(markers_number, by=c("ct1", "ct2", "level") ) %>%

		# # Filter markers
		# filter(rank < `n markers`) %>%
		#
		# # Filter markers in upper levels that have been selected for lower levels
		# inner_join(
		# 	(.) %>%
		# 		distinct(symbol, level) %>%
		# 		group_by(symbol) %>%
		# 		arrange(level %>% desc) %>%
		# 		slice(1) %>%
	# 		ungroup
	# ) %>%

	# Write table
	{
		(.)  %>%	write_csv(sprintf("%s/markers_pass%s.csv", out_dir, pass))
		(.)
	}
}

create_ref = function(ref, markers){

	sample_blacklist = c("666CRI", "972UYG", "344KCP", "555QVG", "370KKZ", "511TST", "13816.11933", "13819.11936", "13817.11934", "13818.11935", "096DQV", "711SNV")


	ref %>%
		inner_join(
			(.) %>%
				distinct(sample) %>%
				filter( !grepl(sample_blacklist %>% paste(collapse="|"), sample))
		) %>%

		# Get house keeping and markwrs
		left_join(markers %>% distinct(level, symbol, ct1, ct2, rank) %>% filter(rank < 500)) %>%
		filter(`house keeping` | rank %>% is.na %>% `!`) %>%

		# Filter out symbol if present both in markers and house keeping (strange but happens)
		anti_join(
			(.) %>% filter(`house keeping` & (ct1 %>% is.na %>% `!`)) %>% distinct(symbol)
		) %>%

		select(  -contains("idx")) %>%
		mutate(`read count` = `read count` %>% as.integer)
}

get_input_data = function(markers, reps, pass){

	sample_blacklist = c("666CRI", "972UYG", "344KCP", "555QVG", "370KKZ", "511TST", "13816.11933", "13819.11936", "13817.11934", "13818.11935", "096DQV", "711SNV")


	list(
		reference =

			ref %>%
			inner_join(
				(.) %>%
					distinct(sample) %>%
					filter( !grepl(sample_blacklist %>% paste(collapse="|"), sample))
			) %>%

			# Filter FANTOM
			#filter(`Data base` != "FANTOM5") %>%

		# inner_join(
		# 	(.) %>%
		# 		distinct(sample, `Cell type category`) %>%
		# 		group_by( `Cell type category`) %>%
		# 		slice(1) %>%
		# 		ungroup
		# ) %>%

		# # Add extended house keeping genes to match the number of markers (hopefully this will solve some exposure divergencies when I have many markers. I have tried with weighting the sampling but did not work)
		# left_join(
		# 	(.) %>%
		# 		filter(level ==2) %>%
		# 		distinct(`Cell type category`, symbol, `house keeping`, lambda, sigma_raw) %>%
		# 		group_by(symbol, `house keeping`) %>%
		# 		summarise(sd = lambda %>% sd, lambda_avg = lambda %>% mean, sigma_raw_sd = sigma_raw %>% sd, sigma_raw_avg = sigma_raw %>% mean) %>%
		# 		ungroup() %>%
		# 		arrange(sd, lambda_avg %>% desc) %>%
		# 		mutate(`house keeping extended rank` = 1:n()) %>%
		# 		slice(1:(markers %>% distinct(level, symbol) %>% nrow))
		# ) %>%
		# 	mutate(`house keeping` = `house keeping extended rank` %>% is.na %>% `!`) %>%
		# 	mutate(lambda = ifelse(`house keeping extended rank` %>% is.na %>% `!`, lambda_avg, lambda)) %>%
		# 	mutate(sigma_raw = ifelse(`house keeping extended rank` %>% is.na %>% `!`, sigma_raw_avg, sigma_raw)) %>%
		# 	select(-sd,-lambda_avg,- sigma_raw_sd,-sigma_raw_avg,- `house keeping extended rank`) %>%

		# Get house keeping and markwrs
		left_join(markers %>% distinct(level, symbol, ct1, ct2, rank, `n markers`) %>% filter(rank < 500)) %>%
			filter(`house keeping` | rank %>% is.na %>% `!`) %>%

			# Filter out symbol if present both in markers and house keeping (strange but happens)
			anti_join(
				(.) %>% filter(`house keeping` & (ct1 %>% is.na %>% `!`)) %>% distinct(symbol)
			) %>%

			# {
			# 	bind_rows(
			# 		(.) %>% inner_join( markers %>% distinct(level, symbol, ct1, ct2)),
			# 		(.) %>% filter(`house keeping`)
			# 	)
			# } %>%

			# decrease number of house keeping
			# anti_join(
			# 	(.) %>% filter(`house keeping`) %>% distinct(symbol) %>% slice(500) #sample_frac(0.7)
		# ) %>%

		select(  -contains("idx")) %>%
			mutate(`read count` = `read count` %>% as.integer)	,
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
			save(input, file=sprintf("dev/input_test_pass_%s.RData", pass))
			(.)
		}


}

# mix_base = readRDS("dev/mix_base.RDS")
# my_ref = 	mix_base %>% distinct(sample, `Cell type category`, level)
#
# set.seed(123)
#
# get_mix_source = function(){
# 	{
#
# 		# This function goes thought nodes and grubs names of the cluster
# 		gn = function(node) {
# 			if(length(node$children) > 0) {
#
# 				result =
#
# 					# Get cildren names
# 					#tibble(parent = node$name, children = foreach(cc = node$children, .combine = c) %do% {cc$name}) %>%
# 					foreach(cc = node$children, .combine = c) %do% {cc$name} %>%
#
# 					#create cmbinations
# 					combn(m = 2) %>%
# 					t %>%
# 					as_tibble %>%
# 					mutate(parent = node$name, level = node$level )
#
# 				# Merge results with other nodes
# 				result %<>%
# 					bind_rows(
# 						foreach(nc = node$children, .combine = bind_rows) %do% {gn(nc)}
# 					)
#
# 				return (result)
# 			}
# 		}
#
# 		gn(Clone(ARMET::tree))
# 	} %>%
# 		mutate(`#` = 1:n()) %>%
#
# 		# Make more runs
# 		right_join(
# 			tibble(`#` = (1: ((.) %>% nrow)) %>% rep(reps))  %>%
# 				mutate(run = 1:n())
# 		) %>%
# 		select(-`#`) %>%
#
# 		#multidplyr::partition(run) %>%
# 		group_by(run) %>%
# 		do({
# 			`%>%` = magrittr::`%>%`
#
# 			cc = (.)
#
# 			bind_rows(
# 				my_ref %>%
# 					filter(`Cell type category` == (cc %>% pull(V1)) & level == (cc %>% pull(level))) %>%
# 					sample_n(1) %>%
# 					distinct(sample, `Cell type category`),
# 				my_ref %>%
# 					filter(`Cell type category` == (cc %>% pull(V2))  & level == (cc %>% pull(level))) %>%
# 					sample_n(1) %>%
# 					distinct(sample, `Cell type category`)
# 			) %>%
# 				mutate(run = cc %>% distinct(run) %>% pull(1)) %>%
# 				mutate(level = cc %>% distinct(level) %>% pull(1))
# 		}) %>%
# 		#multidplyr::collect() %>%
# 		ungroup() %>%
#
# 		# Again solving the problem with house keeping genes
# 		left_join(mix_base  %>% distinct(`symbol`, `read count normalised bayes`, `Cell type category`, sample, level, `house keeping`))  %>%
#
# 		# Add mix_sample
# 		left_join(
# 			(.) %>%
# 				group_by(run) %>%
# 				do(
# 					(.) %>%
# 						distinct(`symbol`, `read count normalised bayes`, `Cell type category`, run, level) %>%
# 						spread(`Cell type category`, `read count normalised bayes`) %>%
# 						drop_na %>%
# 						mutate(pair = names((.))[4:5] %>% paste(collapse=" ")) %>%
# 						setNames(c("symbol", "run", "level", "1", "2", "pair")) %>%
# 						mutate( `read count mix` = ( (`1` + `2`) / 2 ) %>% as.integer ) %>%
# 						unite(sample_mix, c("run", "pair"), remove = F)
# 				) %>%
# 				ungroup() %>%
# 				distinct(symbol, run, level, pair, sample_mix,  `read count mix`)
# 		) %>%
#
# 		# Eliminate duplicated of difference levels for example house keepng genes
# 		group_by(run, symbol, sample) %>%
# 		arrange(level) %>%
# 		slice(1) %>%
# 		ungroup
# }
#
# mix_source = get_mix_source()
#
# saveRDS(mix_source, file="dev/mix_source.rds")

# (xx %>%
# 		#filter(symbol %in% (read_csv("docs/hk_600.txt", col_names = FALSE) %>% pull(1))) %>%
# 		#inner_join((.) %>% distinct(symbol) %>% head(n=300)) %>%
# 		aggregate_duplicated_gene_symbols(value_column = "read count normalised bayes") %>%
# 		add_MDS_components(feature_column  = "symbol", elements_column  = "sample", value_column = "read count normalised bayes", components_list = list(1:2, 3:4, 5:6)) %>%
# 		#rotate_MDS_components(rotation_degrees = -50) %>%
# 		select(contains("Dimension"), sample,  `Cell type formatted`, `Data base`) %>%
# 		distinct() %>%
# 		GGally::ggpairs(columns = 1:6, ggplot2::aes(colour=`Data base`, label=`Cell type formatted`))
# 	) %>% plotly::ggplotly()


# foreach(n_markers = c(1:10, seq(15, 50, 5), seq(50, 100, 10))) %do% {
#

s = sample(size = 1, 1:99999999)

file_name = sprintf("%s/markers_%s_%s.RData", out_dir, n_markers, s)
#
# 	if(file.exists(file_name) %>% `!`) {
# 		try(

set.seed(s)

ARMET_tc(
	readRDS("dev/mix_source.rds") %>%
		distinct(symbol, sample_mix, `read count mix`) %>%
		drop_na %>%
		spread(`sample_mix`, `read count mix`) %>%
		drop_na %>%
		gather(sample_mix, `read count mix`, -symbol) %>%
		spread(`symbol`, `read count mix`) %>%
		rename(sample = sample_mix) %>%
		inner_join((.) %>% distinct(sample) %>% sample_n(50)),
	iterations = 250,
	n_markers = ARMET::ARMET_ref %>% distinct(ct1, ct2) %>% mutate(`n markers` = n_markers),
	full_bayes = T,
	cores = 10, levels = levels
) %>%
	saveRDS(file=file_name)
# 		)
# 	}
# }



#qsub -l nodes=1:ppn=12,mem=32gb,walltime=18:00:00 dev/job_torque.sh -F "20 fist_run"

# # Plot results
# res_dir = "dev/feature_selection_fist_run_normal_queue/"
#
# res =
# 	dir(res_dir, full.names = T) %>%
# 		map_dfr(
# 			~ .x %>%
# 				readRDS() %$%
# 				proportions %>%
# 				mutate(n_markers = .x)
# 		) %>%
# 	separate(sample, c("run", "pair"), sep="_", extra="merge") %>%
# 	separate(pair, c("ct1", "ct2"), sep=" ") %>%
# 	filter(`Cell type category` == ct1 | `Cell type category` == ct2) %>%
# 	separate(n_markers, c("path", "n_markers"), sep="markers_") %>%
# 	separate(n_markers, c("n_markers", "extension"), sep="\\.") %>%
# 	mutate(CI = .upper - .lower) %>%
# 	mutate(error = .value %>% `+` (-0.5) %>% abs)
#
# res %>%
# 	filter(level ==4) %>%
# 	mutate(n_markers = n_markers %>% as.integer) %>%
# 	ggplot(aes(x=n_markers, y=error, color=`Cell type category`)) +
# 	geom_point() +
# 	geom_smooth(method = "lm")
#
# (res %>%
# 	filter(level ==2) %>%
# 	mutate(n_markers = n_markers %>% as.integer) %>%
# 	ggplot(aes(x=n_markers, y=error, color=`Cell type category`, size=CI, Q=Q, C=C)) +
# 	geom_point() +
# 	geom_smooth(method = "lm")
# ) %>% plotly::ggplotly()

# Create reference data set


# get_markers_number(5, res_4, n_markers_4) %>%
# 	get_markers_df(5) %>%
# 	get_input_data(reps = reps, pass = 5) %$%
# 	reference %>%
# 	write_csv(sprintf("%s/reference.csv", out_dir))
#
# # Statistics
#
# (
# 	list(
# 		n_markers_0 %>% mutate(pass = 0),
# 		n_markers_1 %>% mutate(pass = 1),
# 		n_markers_2 %>% mutate(pass = 2),
# 		n_markers_3 %>% mutate(pass = 3),
# 		n_markers_4 %>% mutate(pass = 4)
# 	) %>%
# 		bind_rows() %>%
# 		ggplot(aes(x=pass, y=`n markers`, color=pair)) + geom_point() + geom_line() + theme(legend.title = element_blank())
# ) %>%
# 	plotly::ggplotly()
#
# (
# 	list(
# 		n_markers_0 %>% mutate(pass = 0),
# 		n_markers_1 %>% mutate(pass = 1),
# 		n_markers_2 %>% mutate(pass = 2),
# 		n_markers_3 %>% mutate(pass = 3),
# 		n_markers_4 %>% mutate(pass = 4)
# 	) %>%
# 		bind_rows() %>%
# 		ggplot(aes(x=pass, y=`error mean relative mean`, color=pair)) + geom_point() + geom_line() + theme(legend.title = element_blank())
# ) %>%
# 	plotly::ggplotly()
