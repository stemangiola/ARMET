# Create mix_base within NA - imputed values

library(tidyverse)
library(magrittr)
library(purrr)
library(furrr)
library(data.tree)
library(foreach)
library(ARMET)

n_cores = 20
S = 30

my_theme =
	theme_bw() +
	theme(
		panel.border = element_blank(),
		axis.line = element_line(),
		panel.grid.major = element_line(size = 0.2),
		panel.grid.minor = element_line(size = 0.1),
		text = element_text(size = 12),
		legend.position = "bottom",
		aspect.ratio = 1,
		strip.background = element_blank(),
		axis.title.x  = element_text(margin = margin(
			t = 10,
			r = 10,
			b = 10,
			l = 10
		)),
		axis.title.y  = element_text(margin = margin(
			t = 10,
			r = 10,
			b = 10,
			l = 10
		)),
		axis.text.x = element_text(angle = 90, hjust = 1)
	)

noise_test = function(which_is_up_down, do_regression = F) {
	# intercept = rnorm(16)
	intercept = rep(1, 16)
	intercept = intercept - sum(intercept) / length(intercept)
	slope = rep(0, 16)
	slope[which_is_up_down[1]] = 2
	slope[which_is_up_down[2]] = -4
	slope[which_is_up_down[3]] = 1
	slope[which_is_up_down[4]] = 0.5
	alpha = matrix(intercept %>%	c(slope), ncol = 2)

	mix = mix_base %>% generate_mixture(15, alpha)

	rr =
		mix %>%
		select(run, symbol, `count mix`, covariate_2) %>%
		mutate(`count mix` = as.integer(`count mix`), run = as.character(run)) %>%
		spread(symbol, `count mix`) %>%
		rename(sample = run) %>%
		ARMET_tc(
			~ covariate_2,
			do_regression = do_regression,
			iterations = ifelse(do_regression, 2400, 700),
			sampling_iterations = ifelse(do_regression, 2000, 300)
		)

	list(mix = mix, result = rr)
	#%>% saveRDS(sprintf("dev/test_student_noisless_%s", which_is_up_down %>% paste(collapse="_")))
}

mix_base = readRDS("dev/mix_base.RDS") %>% filter(level==3)

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
				gather(`Cell type category`, alpha, -run)
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

		left_join(cell_type_proportions %>% distinct(run, covariate_2)) %>%

		# Add proportions
		add_attr(cell_type_proportions, "proportions")

}

which_is_up_down = 1:16 %>% map( ~ c(
	.x,
	(.x + 4) %>% ifelse(. > 16, . - 16, .),
	(.x + 8) %>% ifelse(. > 16, . - 16, .),
	(.x + 12) %>% ifelse(. > 16, . - 16, .)
))

# which_is_up_down %>%
# 	map( ~ .x %>% noise_test) %>%
# 	saveRDS("dev/test_dirichlet_with_noise.rds")

which_is_up_down %>%
	map( ~ .x %>% noise_test(do_regression = T)) %>%
	saveRDS("dev/test_dirichlet_with_noise_regression.rds")
