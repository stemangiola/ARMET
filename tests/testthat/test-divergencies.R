context('Test if divergent')
library(rlang)
library(tidyverse)
library(ARMET)
library(data.tree)
library(foreach)
library(magrittr)

mix_base = readRDS("~/PhD/deconvolution/ARMET/dev/mix_base.RDS")
my_ref = 	mix_base %>% distinct(sample, `Cell type category`, level)
reps = 10
set.seed(123)

# %>%
# 	filter( ! grepl(sample_blacklist %>% paste(collapse="|"), sample))

#
marker_df  = ARMET::ARMET_ref

get_mix_source = function(){
	{

		# This function goes thought nodes and grubs names of the cluster
		gn = function(node) {
			if(length(node$children) > 0) {

				result =

					# Get cildren names
					#tibble(parent = node$name, children = foreach(cc = node$children, .combine = c) %do% {cc$name}) %>%
					foreach(cc = node$children, .combine = c) %do% {cc$name} %>%

					#create cmbinations
					combn(m = 2) %>%
					t %>%
					as_tibble %>%
					mutate(parent = node$name, level = node$level )

				# Merge results with other nodes
				result %<>%
					bind_rows(
						foreach(nc = node$children, .combine = bind_rows) %do% {gn(nc)}
					)

				return (result)
			}
		}

		gn(Clone(ARMET::tree))
	} %>%
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
					filter(`Cell type category` == (cc %>% pull(V1)) & level == (cc %>% pull(level))) %>%
					sample_n(1) %>%
					distinct(sample, `Cell type category`),
				my_ref %>%
					filter(`Cell type category` == (cc %>% pull(V2))  & level == (cc %>% pull(level))) %>%
					sample_n(1) %>%
					distinct(sample, `Cell type category`)
			) %>%
				mutate(run = cc %>% distinct(run) %>% pull(1)) %>%
				mutate(level = cc %>% distinct(level) %>% pull(1))
		}) %>%
		#multidplyr::collect() %>%
		ungroup() %>%

		# Again solving the problem with house keeping genes
		left_join(mix_base  %>% distinct(`symbol`, `read count normalised bayes`, `Cell type category`, sample, level, `house keeping`))  %>%

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
}

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

mix_source = get_mix_source()

# infer
res =
  ARMET_tc(
    mix_source %>%
      inner_join( (.) %>% distinct(sample_mix) %>% slice(1:2) ) %>%
      distinct(symbol, sample_mix, `read count mix`) %>%
      drop_na %>%
      spread(`sample_mix`, `read count mix`) %>%
      drop_na %>%
      gather(sample_mix, `read count mix`, -symbol) %>%
      spread(`symbol`, `read count mix`) %>%
      rename(sample = sample_mix),
    n_markers = get_markers_number(0,min_n_samples = 10, NULL, NULL) %>% select(-level) %>% drop_na,
    full_bayes = T,
    cores = 5,
    iterations = 250 #,sampling_iterations = 2
  )

test_that("Test data frame",{ expect_equal( ncol(ttSc::counts_sc), 6 ) })

test_that("Create tt object from tibble",{

  expect_equal( ncol(tt), 11 )

  expect_equal( typeof(attr(tt, "parameters")), "list")

})
