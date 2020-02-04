# Create mix_base withn NA - imputed values

library(tidyverse)
library(magrittr)
library(purrr)
library(furrr)
library(data.tree)
library(foreach)
library(ARMET)
library(gamlss)
source("~/PhD/deconvolution/ARMET/R/utils.R")
source("~/PostDoc/ppcSeq/R/do_parallel.R")
n_cores = 20
S = 30

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
		strip.background = element_blank(),
		axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
		axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
		axis.text.x = element_text(angle = 90, hjust = 1)
	)

noiseles_test = function(which_is_up_down){
	# intercept = rnorm(16)
	intercept = rep(1, 16)
	intercept = intercept - sum(intercept) / length(intercept)
	slope = rep(0, 16)
	slope[which_is_up_down[1]] = 2
	slope[which_is_up_down[2]] = -2
	alpha = matrix(intercept %>%	c(slope), ncol = 2)

	mix = mix_base %>% generate_mixture(15, alpha)

	rr =
		mix %>%
		select(run, symbol, `count mix`, covariate_2) %>%
		mutate(`count mix` = as.integer(`count mix`), run = as.character(run) ) %>%
		spread(symbol, `count mix`) %>%
		rename(sample = run) %>%
		ARMET_tc(~covariate_2, do_regression = T, iterations = 700, sampling_iterations = 200	)

	list(mix = mix, result = rr)
	#%>% saveRDS(sprintf("dev/test_student_noisless_%s", which_is_up_down %>% paste(collapse="_")))
}

mix_base = readRDS("dev/mix_base_noiseless.RDS") %>% filter(level == 3)

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

which_is_up_down = 1:16 %>% map(~ c(.x, (.x + 4) %>% ifelse(. > 16, .- 16, .)))

which_is_up_down %>%
	map(~ .x %>% noiseles_test) %>%
	saveRDS("dev/test_student_noisless.rds")

res = readRDS("dev/test_student_noisless.rds")

# my_res = res[[1]]
#
# my_res$mix %>% attr("proportions") %>%
# 	ggplot(aes(x = covariate_2, y = p, color=factor(alpha_2))) + geom_point() + geom_smooth() + facet_wrap(~`Cell type category`)

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
	ggplot(aes(tpr, x=`Cell type category`)) +
	geom_boxplot(outlier.shape = NA) +
	geom_jitter() +
	my_theme

# Calculate fpr for regression
res %>%
	map_dfr(
		~
			# Integrate
			.x$result$proportions %>%
			filter(level == 3) %>%
			select(`Cell type category`, contains("alpha2")) %>%
			distinct() %>%
			left_join(.x$mix %>% attr("proportions") %>% distinct(`Cell type category`, alpha_2) ) %>%
			drop_na  %>%

			# Calculate
			mutate(fp = alpha_2 == 0 & (.lower_alpha2 * .upper_alpha2 ) > 0) %>%
			mutate(fn = alpha_2 != 0 & (.lower_alpha2 * .upper_alpha2 ) < 0)

	) %>%
	group_by(`Cell type category`) %>%
	summarise(fpr = sum(fp) / n(), fnr = sum(fn) / n()) %>%
	gather(which, rate, c("fpr","fnr")) %>%
	ggplot(aes(y=rate, x=`Cell type category`)) +
	geom_point() +
	facet_wrap(~which) +
	my_theme

# Example of one run
comparison_truth %>%
	ggplot(aes(x = p, y = .value, color=`Cell type category`)) +
	geom_abline(intercept = 0 ,slope=1) +
	geom_smooth(method = "lm") +
	geom_point() +
	geom_errorbar(aes(ymin=.value.lower, ymax=.value.upper))

# Same run with Cibersort

noiseles_test_cibersort = function(which_is_up_down, ref){
	# intercept = rnorm(16)
	intercept = rep(1, 16)
	intercept = intercept - sum(intercept) / length(intercept)
	slope = rep(0, 16)
	slope[which_is_up_down[1]] = 2
	slope[which_is_up_down[2]] = -2
	alpha = matrix(intercept %>%	c(slope), ncol = 2)

	mix =
		mix_base %>% generate_mixture(15, alpha)

	rr =
		mix %>%
		mutate(run = as.character(run)) %>%
		ttBulk::deconvolve_cellularity(run, symbol, `count mix`, reference = ref) %>%
		dplyr::select(run, covariate_2, contains("type")) %>%
		distinct() %>%
		pivot_longer(cols = contains("type"), names_to = "Cell type category", values_to = ".value",  names_prefix = "type: ") %>%

		# perform test
		nest(data = -c(`Cell type category`)) %>%
		mutate(CI = data %>% map(~ .x %>% gamlss::gamlss(.value ~ covariate_2,  family = BEZI, data = ., trace = F) %>%
											summary %>% `[` (2, 1:2) %>% as_tibble(rownames="rn") %>% spread(rn, value) %>%
											rename(	alpha_2 = Estimate) %>%
											mutate(
												.lower_alpha2 = alpha_2 - `Std. Error`,
												.upper_alpha2 = alpha_2 + `Std. Error`
											) %>%
											dplyr::select(-`Std. Error`))
				) %>%
		unnest(cols = c(data, CI))

	list(mix = mix, result = rr)
	#%>% saveRDS(sprintf("dev/test_student_noisless_%s", which_is_up_down %>% paste(collapse="_")))
}

# get reference from ARMET
ref =
	res[[1]]$result$signatures[[3]] %>%
	filter(!`Cell type category` %in% c("house_keeping", "query")) %>%
	distinct(`Cell type category`, symbol, lambda_log) %>%
	mutate(count = exp(lambda_log)) %>%
	dplyr::select(-lambda_log) %>%
	spread(`Cell type category`, count) %>%
	as_matrix(rownames="symbol")


res_cibersort =
	which_is_up_down %>%
	map(~ .x %>% noiseles_test_cibersort(ref))

res_cibersort[[1]]$mix %>% attr("proportions") %>%
	mutate(run = run %>% as.character) %>%
	dplyr::select(-contains("alpha")) %>%
	left_join(res_cibersort[[1]]$result) %>%
	ggplot(aes(x = p, y = .value, color=`Cell type category`)) +
	geom_abline(intercept = 0 ,slope=1) +
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
	ggplot(aes(tpr, x=`Cell type category`)) +
	geom_boxplot(outlier.shape = NA) +
	geom_jitter() +
	my_theme

# Calculate fpr for regression
res_cibersort %>%
	map_dfr(
		~
			# Integrate
			.x$result %>%
			dplyr::select(`Cell type category`, contains("alpha2")) %>%
			distinct() %>%
			left_join(.x$mix %>% attr("proportions") %>% distinct(`Cell type category`, alpha_2) ) %>%
			drop_na  %>%

			# Calculate
			mutate(fp = alpha_2 == 0 & (.lower_alpha2 * .upper_alpha2 ) > 0) %>%
			mutate(fn = alpha_2 != 0 & (.lower_alpha2 * .upper_alpha2 ) < 0)

	) %>%
	group_by(`Cell type category`) %>%
	summarise(fpr = sum(fp) / n(), fnr = sum(fn) / n()) %>%
	gather(which, rate, c("fpr","fnr")) %>%
	ggplot(aes(y=rate, x=`Cell type category`)) +
	geom_point() +
	facet_wrap(~which) +
	my_theme
