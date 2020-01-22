# Create mix_base withn NA - imputed values

library(tidyverse)
library(magrittr)
library(purrr)
library(furrr)
library(data.tree)
library(foreach)
library(ARMET)
source("~/PhD/deconvolution/ARMET/R/utils.R")
source("~/PostDoc/ppcSeq/R/do_parallel.R")
n_cores = 20
S = 30



mix_base = readRDS("dev/mix_base.RDS") %>% filter(level == 3)

intercept = c(	-0.9516997,	0.3748425,	0.2348878,	-0.4229247,	-0.9586268,	-3.1402425,	0.3323535,-0.7435126,	2.9877128,	-1.2978599,	4.3127303,	0.5863942,	0.3885436,	1.3248586, 0.3885436,	1.3248586)
intercept = intercept - sum(intercept) / length(intercept)
alpha = matrix(intercept %>%	 	c(0 , 0 , 0, 2, 0, 2, 0, 0, 2, 0, 0, 2, 0, 2, 0, 0), ncol = 2)

generate_mixture = function(.data, samples_per_condition, alpha) {
	add_attr = function(var, attribute, name) {
		attr(var, name) <- attribute
		var
	}

	logsumexp <- function (x) {
		y = max(x)
		y + log(sum(exp(x - y)))
	}

	softmax <- function (x) {
		exp(x - logsumexp(x))
	}

	S = samples_per_condition * 2

	X = matrix(rep(1, S) %>%
						 	c(rep(0, S / 2)) %>%
						 	c(rep(1, S / 2)),
						 ncol = 2)

	samples_per_run =
		foreach(r = 1:S, .combine = bind_rows) %do% {
			.data %>%
				distinct(`Cell type category`, sample) %>%
				group_by(`Cell type category`) %>%
				sample_n(1) %>%
				ungroup() %>%
				mutate(run = r)
		}

	ct_names = .data %>% distinct(`Cell type category`) %>% pull(1)

	alpha_df = alpha %>% as.data.frame %>% setNames(sprintf("alpha_%s", 1:2)) %>% mutate(`Cell type category`  = ct_names)

	cell_type_proportions =
		# Choose samples
		samples_per_run %>%

		# Choose proportions
		left_join(
			# Decide theoretical, noise-less proportions for each sample
			X %*% t(alpha) %>%
				apply(1, softmax) %>%
				t %>%
				`*` (40) %>%
				as.data.frame() %>%
				as_tibble() %>%
				setNames(ct_names) %>%
				mutate(run = 1:n()) %>%
				gather(`Cell type category`, alpha,-run)
		) %>%

		# Add X
		left_join(X %>% as.data.frame %>% setNames(sprintf("covariate_%s", 1:2)) %>% mutate(run = 1:n())) %>%

		# Add alpha
		left_join(alpha_df) %>%

		group_by(run) %>%
		mutate(p = gtools::rdirichlet(1, alpha)) %>%
		ungroup()

	# Add counts
	dirichlet_source =
		cell_type_proportions %>%
		left_join(.data, by = c("Cell type category", "sample"))

	# Make mix
	dirichlet_source %>%
		mutate(c = `count normalised bayes` * p) %>%
		group_by(run, symbol) %>%
		summarise(`count mix` = c %>% sum) %>%
		ungroup %>%

		left_join(cell_type_proportions %>% distinct(run, covariate_2 )	) %>%

		# Add proportions
		add_attr(cell_type_proportions, "proportions")

}

mix = mix_base %>% generate_mixture(15, alpha)

rr2 =
	mix %>%
	select(run, symbol, `count mix`, covariate_2) %>%
	mutate(`count mix` = as.integer(`count mix`), run = as.character(run) ) %>%
	spread(symbol, `count mix`) %>%
	rename(sample = run) %>%
	ARMET_tc(~covariate_2, do_regression = T	)

mix %>% attr("proportions") %>%
	ggplot(aes(x = covariate_2, y = p, color=factor(alpha_2))) + geom_point() + geom_smooth() + facet_wrap(~`Cell type category`)

# See if result match
rr2$proportions %>%
	filter(level == 3) %>%
	select(`Cell type category`, contains("value_alpha")) %>%
	distinct() %>%
	left_join(mix %>% attr("proportions") %>% distinct(`Cell type category`, alpha_2) ) %>%
	drop_na
