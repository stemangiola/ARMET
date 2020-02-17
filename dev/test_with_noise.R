# Create mix_base withn NA - imputed values

library(tidyverse)
library(magrittr)
library(purrr)
library(furrr)
library(data.tree)
library(foreach)
library(ARMET)
library(broom)
# source("~/PostDoc/ppcSeq/R/do_parallel.R")
# n_cores = 20
# S = 30

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

which_is_up_down %>%
	map( ~ .x %>% noise_test) %>%
	saveRDS("dev/test_dirichlet_with_noise.rds")

which_is_up_down %>%
	map( ~ .x %>% noise_test(do_regression = T)) %>%
	saveRDS("dev/test_dirichlet_with_noise_regression.rds")

res = readRDS("dev/test_dirichlet_with_noise.rds")

res_regression = readRDS("dev/test_dirichlet_with_noise_regression.rds")

# my_res = res[[1]]
#
# my_res$mix %>% attr("proportions") %>%
# 	ggplot(aes(x = covariate_2, y = p, color=factor(alpha_2))) + geom_point() + geom_smooth() + facet_wrap(~`Cell type category`)

# Example of one run
res[[1]]$mix %>% attr("proportions") %>%
	mutate(sample = run %>% as.character) %>%
	dplyr::select(-contains("alpha")) %>%
	left_join(res[[1]]$result$proportions) %>%
	ggplot(aes(x = p, y = .value, color = `Cell type category`)) +
	geom_abline(intercept = 0 , slope = 1) +
	geom_smooth(method = "lm") +
	geom_errorbar(aes(ymin = .value.lower, ymax = .value.upper), alpha = 0.2) +
	geom_point()

res_regression[[8]]$result$fit[[3]] %>% tidybayes::gather_draws(alpha_3[A,C]) %>% filter(A==2) %>% ggplot(aes(.value, color=factor(C))) + geom_density() + geom_vline(xintercept = 0.232)

# Calculate tpr
res %>%
	map_dfr(
		~
			# Integrate
			.x$result$proportions %>%
			filter(level == 3) %>%
			rename(run = sample) %>%
			select(run, `Cell type category`, .value, .value.lower, .value.upper) %>%
			left_join(
				.x$mix %>% attr("proportions") %>% mutate(run = run %>% as.character)
			) %>%

			# Calculate
			mutate(inside = p > .value.lower & p < .value.upper) %>%
			count(`Cell type category`, inside) %>%
			spread(inside, n) %>%
			replace(., is.na(.), 0) %>%
			mutate(tpr = `TRUE` / (`TRUE` + `FALSE`))
	) %>%
	ggplot(aes(tpr, x = `Cell type category`)) +
	geom_boxplot(outlier.shape = NA) +
	geom_jitter() +
	my_theme

# Calculate fpr for regression
regression_ARMET =
	res_regression %>%
	map_dfr(
		~
			# Integrate
			.x$result$proportions %>%
			filter(level == 3) %>%
			select(`Cell type category`, contains("alpha2"), zero) %>%
			distinct() %>%
			left_join(
				.x$mix %>% attr("proportions") %>% distinct(`Cell type category`, alpha_2)
			) %>%
			drop_na  %>%

			# Subtract relative zero
			mutate(
				.value_alpha2 = .value_alpha2 - zero,
				.lower_alpha2 = .lower_alpha2 - zero,
				.upper_alpha2 = .upper_alpha2 - zero
			) %>%

			# Calculate
			mutate(fp = alpha_2 == 0 &
						 	(.lower_alpha2 * .upper_alpha2) > 0) %>%
			mutate(fn = alpha_2 != 0 & (.lower_alpha2 * .upper_alpha2) < 0)

	)

regression_ARMET %>%
	group_by(`Cell type category`, alpha_2) %>%
	summarise(fpr = sum(fp) / n(), fnr = sum(fn) / n()) %>%
	gather(which, rate, c("fpr", "fnr")) %>%
	ggplot(aes(y = rate, x = `Cell type category`, color = factor(alpha_2))) +
	geom_jitter() +
	facet_wrap( ~ which) +
	my_theme


# get reference from ARMET
ref =
	res[[1]]$result$signatures[[3]] %>%
	filter(!`Cell type category` %in% c("house_keeping", "query")) %>%
	distinct(`Cell type category`, symbol, lambda_log) %>%
	mutate(count = exp(lambda_log)) %>%
	dplyr::select(-lambda_log) %>%
	spread(`Cell type category`, count) %>%
	ttBulk::as_matrix(rownames = "symbol")

noise_test_ttBulk = function(which_is_up_down, ref, method = "cibersort") {
	# intercept = rnorm(16)
	intercept = rep(1, 16)
	intercept = intercept - sum(intercept) / length(intercept)
	slope = rep(0, 16)
	slope[which_is_up_down[1]] = 2
	slope[which_is_up_down[2]] = -4
	slope[which_is_up_down[3]] = 1
	slope[which_is_up_down[4]] = 0.5
	alpha = matrix(intercept %>%	c(slope), ncol = 2)

	mix =
		mix_base %>% generate_mixture(15, alpha)

	rr =
		mix %>%
		mutate(run = as.character(run)) %>%
		ttBulk::deconvolve_cellularity(run, symbol, `count mix`, reference = ref, method = method) %>%
		dplyr::select(run, covariate_2, contains(method)) %>%
		distinct() %>%
		pivot_longer(
			cols = contains(method),
			names_to = "Cell type category",
			values_to = ".value",
			names_prefix = sprintf("%s: ", method)
		) %>%

		# perform test
		nest(data = -c(`Cell type category`)) %>%
		mutate(
			CI = data %>% map(
				~ .x %>% gamlss::gamlss(
					.value ~ covariate_2,
					family = gamlss.dist::BEZI,
					data = .,
					trace = F
				) %>%
					summary %>% `[` (2, c(1, 2, 4)) %>% as_tibble(rownames = "rn") %>% spread(rn, value) %>%
					rename(alpha_2 = Estimate) %>%
					mutate(
						.lower_alpha2 = alpha_2 - `Std. Error`,
						.upper_alpha2 = alpha_2 + `Std. Error`
					) %>%
					dplyr::select(-`Std. Error`)
			)
		) %>%
		unnest(cols = c(data, CI))

	list(mix = mix, result = rr)
	#%>% saveRDS(sprintf("dev/test_student_noisless_%s", which_is_up_down %>% paste(collapse="_")))
}


# Same run with CIBERSORT

res_cibersort =
	which_is_up_down %>%
	map( ~ .x %>% noise_test_ttBulk(ref))

res_cibersort[[1]]$mix %>% attr("proportions") %>%
	mutate(run = run %>% as.character) %>%
	dplyr::select(-contains("alpha")) %>%
	left_join(res_cibersort[[1]]$result) %>%
	ggplot(aes(x = p, y = .value, color = `Cell type category`)) +
	geom_abline(intercept = 0 , slope = 1) +
	geom_smooth(method = "lm") +
	geom_point()

# Calculate tpr
res_cibersort %>%
	map_dfr(
		~
			# Integrate
			.x$result %>%
			rename(run = sample) %>%
			select(run, `Cell type category`, .value) %>%
			left_join(
				.x$mix %>% attr("proportions") %>% mutate(run = run %>% as.character)
			) %>%

			# Calculate
			mutate(inside = p > .value.lower & p < .value.upper) %>%
			count(`Cell type category`, inside) %>%
			spread(inside, n) %>%
			replace(., is.na(.), 0) %>%
			mutate(tpr = `TRUE` / (`TRUE` + `FALSE`))
	) %>%
	ggplot(aes(tpr, x = `Cell type category`)) +
	geom_boxplot(outlier.shape = NA) +
	geom_jitter() +
	my_theme

# Calculate fpr for regression
regression_cibersort =
	res_cibersort %>%
	map_dfr(
		~
			# Integrate
			.x$result %>%
			dplyr::select(`Cell type category`, contains("alpha2"), 	`Pr(>|t|)`) %>%
			distinct() %>%
			left_join(
				.x$mix %>% attr("proportions") %>% distinct(`Cell type category`, alpha_2),
				by = "Cell type category"
			) %>%
			drop_na  %>%

			# Calculate
			mutate(fp = alpha_2 == 0 &
						 	`Pr(>|t|)` < 0.005) %>%
			mutate(fn = alpha_2 != 0 & `Pr(>|t|)` > 0.005)

	)

regression_cibersort %>%
	group_by(`Cell type category`, alpha_2) %>%
	summarise(fpr = sum(fp) / n(), fnr = sum(fn) / n()) %>%
	gather(which, rate, c("fpr", "fnr"))%>%
	ggplot(aes(y = rate, x = `Cell type category`, color = factor( alpha_2))) +
	geom_jitter() +
	facet_wrap( ~ which) +
	my_theme

# Same run with LLSR
res_llsr =
	which_is_up_down %>%
	map( ~ .x %>% noise_test_ttBulk(ref, method="llsr"))

# Calculate fpr for regression
regression_llsr =
	res_llsr %>%
	map_dfr(
		~
			# Integrate
			.x$result %>%
			dplyr::select(`Cell type category`, contains("alpha2"), 	`Pr(>|t|)`) %>%
			distinct() %>%
			left_join(
				.x$mix %>% attr("proportions") %>% distinct(`Cell type category`, alpha_2),
				by = "Cell type category"
			) %>%
			drop_na  %>%

			# Calculate
			mutate(fp = alpha_2 == 0 &
						 	`Pr(>|t|)` < 0.005) %>%
			mutate(fn = alpha_2 != 0 & `Pr(>|t|)` > 0.005)

	)

regression_llsr %>%
	group_by(`Cell type category`, alpha_2) %>%
	summarise(fpr = sum(fp) / n(), fnr = sum(fn) / n()) %>%
	gather(which, rate, c("fpr", "fnr")) %>%
	ggplot(aes(y = rate, x = `Cell type category`, color = factor(alpha_2))) +
	geom_jitter() +
	facet_wrap( ~ which) +
	my_theme

# Create boxplot of errors
all_results =
	res_llsr %>%
		map_dfr(
			~ .x$mix %>% attr("proportions") %>%
				mutate(run = run %>% as.character) %>%
				dplyr::select(-contains("alpha")) %>%
				left_join(.x$result, by = c("Cell type category", "run", "covariate_2"))
		) %>%
		mutate(.value = ifelse(.value == 0, min(.value[.value!=0]), .value)) %>%
		mutate(model = "lm") %>%

		# Cibersort
		bind_rows(
			res_cibersort %>%
			map_dfr(
				~ .x$mix %>% attr("proportions") %>%
					mutate(run = run %>% as.character) %>%
					dplyr::select(-contains("alpha")) %>%
					left_join(.x$result, by = c("Cell type category", "run", "covariate_2"))
			) %>%
				mutate(.value = ifelse(.value == 0, min(.value[.value!=0]), .value)) %>%
				mutate(model = "cibersort")
		) %>%

		# ARMET
		bind_rows(
			res %>%
				map_dfr(
					~ 	.x$mix %>% attr("proportions") %>%
				mutate(sample = run %>% as.character) %>%
				dplyr::select(-contains("alpha")) %>%
				left_join(.x$result$proportions, by = c("Cell type category", "sample"))
				)%>%
				mutate(run = as.character(run)) %>%
				mutate(model = "ARMET")
		) %>%

		# Add error
		mutate(error_logit_space = abs(gtools::logit(p) - gtools::logit(.value))) %>%

		# Add regression
		nest(data = -c(run, `Cell type category`, model)) %>%
		mutate(
			fit = map(data, ~ lm(p ~ .value , data = .x)),
			tidied = map(fit, tidy),
			glanced = map(fit, glance)
		)

# Plot error
all_results %>%
	unnest(data) %>%
	ggplot(aes(x = model, y = error_gtools::logit_space)) +
	geom_boxplot() +
	facet_wrap(~ `Cell type category`)

# Plot fits to linear function R squared
all_results %>%
	unnest(glanced) %>%
	ggplot(aes(x = model, y = adj.r.squared)) +
	geom_boxplot() +
	facet_wrap(~ `Cell type category`)

# Plot fits to linear function standard error
all_results %>%
	unnest(tidied) %>%
	filter(term ==".value") %>%
	ggplot(aes(x = model, y = std.error)) +
	geom_boxplot() +
	facet_wrap(~ `Cell type category`) + scale_y_log10()

# Regression accuracy
regression_ARMET %>% mutate(algorithm="ARMET") %>%
bind_rows(	regression_cibersort %>% mutate(algorithm="cibersort")) %>%
bind_rows(	regression_llsr %>% mutate(algorithm="llsr")) %>%
group_by( alpha_2, algorithm) %>%
summarise(fpr = sum(fp) / n(), fnr = sum(fn) / n()) %>%
gather(which, rate, c("fpr", "fnr")) %>%
ggplot(aes(y = rate, x = algorithm, color = factor(alpha_2))) +
geom_point() +
facet_grid( ~ which) +
my_theme
