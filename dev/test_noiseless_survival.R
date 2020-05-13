# Create mix_base withn NA - imputed values

library(tidyverse)
library(magrittr)
library(purrr)
library(furrr)
library(data.tree)
library(foreach)
library(ARMET)
#source("~/PhD/deconvolution/ARMET/R/utils.R")
#source("~/PostDoc/ppcSeq/R/do_parallel.R")

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

get_alpha = function(slope, which_changing, cell_types){
	
	# Get the alpha matrix
	
	intercept = rep(0, length(cell_types))
	slope_arr = rep(0, length(cell_types))

	slope_arr[which_changing] = slope
	matrix(intercept %>%	c(slope_arr), ncol = 2)
	
}

get_survival_X = function(S){
	readRDS("dev/PFI_all_cancers.rds") %>%
		filter(PFI.2 == 1 & !is.na(PFI.time.2)) %>%
		select(real_days = PFI.time.2 ) %>%
		mutate(real_days = real_days %>% scale(center = F) %>% as.numeric) %>%
		sample_n(S) %>%
		mutate(sample = sprintf("S%s", 1:n())) %>%
		mutate(alive = sample(0:1, n(), replace = T)) %>%
		mutate(days = ifelse(alive==1, real_days/2, real_days) ) %>%
		mutate(intercept = 1)
}

generate_mixture = function(.data, X_df, alpha) {
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
	
	
	X = X_df %>% select(intercept, real_days) %>% nanny::as_matrix()
	
	samples_per_run =
		map_dfr(
			1:nrow(X), ~ 
				.data %>%
				distinct(`Cell type category`, sample) %>%
				group_by(`Cell type category`) %>%
				sample_n(1) %>%
				ungroup() %>%
				mutate(run = .x)
		)
	
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
		left_join(X_df %>% select(-sample) %>% mutate(run = 1:n())) %>%
		
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
		
		left_join(dirichlet_source %>% nanny::subset(run) ) %>%
		
		# Add proportions
		add_attr(cell_type_proportions, "proportions")
	
}

get_noiseless_harmonised = function(){

	mix_base_unharmonized = readRDS("dev/mix_base_noiseless.RDS")

	my_markers =
		ARMET::ARMET_ref %>%

		left_join(ARMET::n_markers, by = c("ct1", "ct2")) %>%
		filter_reference(
			mix_base_unharmonized %>%
				filter(level == 3) %>%
				distinct(`Cell type category`, symbol, `count normalised bayes`) %>%
				spread(symbol, `count normalised bayes`),
			ARMET::n_markers
		) %>% distinct(level, symbol)

	# level 1
	abundance_1 =
		my_markers %>% filter(level == 1) %>%
		left_join(mix_base_unharmonized) %>%
		select(level_2, symbol,  `count normalised bayes 1` =`count normalised bayes`)

	abundance_2 =
		my_markers %>% filter(level == 2) %>%
		left_join(mix_base_unharmonized) %>%
		select(level_3, symbol,  `count normalised bayes 2` =`count normalised bayes`)

	# Now this is noiseless for the ancestor markers so also for ARMET that rely on hierarchy
	mix_base_unharmonized %>%
	filter(level==3) %>%
	left_join(abundance_2) %>%
	left_join(abundance_1) %>%
	mutate(`count normalised bayes 2` = ifelse(`count normalised bayes 1` %>% is.na, `count normalised bayes 2`, `count normalised bayes 1`)) %>%
	mutate(`count normalised bayes` = ifelse(`count normalised bayes 2` %>% is.na, `count normalised bayes`, `count normalised bayes 2`)) %>%
		select(level_2, level_3, level_4, `Cell type category`, level, sample, symbol, `count normalised bayes`, `house keeping`)

}

noiseles_test = function(mix) {
	

	
	rr =
		mix %>%
		mutate(`count mix` = as.integer(`count mix`), run = as.character(run)) %>%
		select(-level) %>%
		ARMET_tc(
			~ censored(days, alive),
			run, symbol, `count mix`,
			prior_survival_time = X_df$real_days %>% as.numeric
		)  %>%
		ARMET_tc_continue(2) %>%
		ARMET_tc_continue(3)
	
	list(mix = mix, result = rr)
	#%>% saveRDS(sprintf("dev/test_student_noisless_%s", which_is_up_down %>% paste(collapse="_")))
}

get_mix = function(mix_base, slope, which_changing, S){
	
	cell_types =  mix_base %>% filter(level ==3) %>% pull(`Cell type category`) %>% unique
	
	alpha = get_alpha(slope, which_changing, cell_types)
	X_df = get_survival_X(S)
	
	mix = mix_base %>% generate_mixture(X_df, alpha)
	
}

mix_base = get_noiseless_harmonised()

which_is_up_down = 1:16

mixes =
	which_is_up_down %>%
	map( ~ get_mix(mix_base, 4, .x, 30))

# mixes %>%
# 	map( ~ noiseles_test(.x)) %>%
# 	saveRDS("dev/test_noisless_survival.rds")
#
# which_is_up_down %>%
# 	map( ~ .x %>% noiseles_test(do_regression = T)) %>%
# 	saveRDS("dev/test_student_noisless_regression.rds")

res = readRDS("dev/test_noisless_survival.rds")

# res = noiseles_test(mix_base, 4, 1, 30)

# res$result %>% test_differential_composition() %>% filter(significant)

# res$result %>% plot_scatter()


# Density
# (res$result$proportions %>%
# 		#filter(level ==3) %>%
# 		unnest(draws) %>%
# 		filter(A == 2) %>%
# 		ggplot(aes(.value, color=`Cell type category`)) +
# 		geom_density() +
# 		facet_wrap(~.variable)
# ) %>% plotly::ggplotly()
# #


# my_res = res[[1]]
#
# my_res$mix %>% attr("proportions") %>%
# 	ggplot(aes(x = covariate_2, y = p, color=factor(alpha_2))) + geom_point() + geom_smooth() + facet_wrap(~`Cell type category`)

# Example of one run
res[[1]]$mix %>% attr("proportions") %>%
	mutate(sample = run %>% as.character) %>%
	dplyr::select(-contains("alpha")) %>%
	left_join(res[[1]]$result$proportions %>% select(-.variable) %>% unnest(proportions)) %>%
	ggplot(aes(x = p, y = .value, color = `Cell type category`)) +
	geom_abline(intercept = 0 , slope = 1) +
	geom_smooth(method = "lm") +
	geom_errorbar(aes(ymin = .value.lower, ymax = .value.upper), alpha = 0.2) +
	geom_point()

# Calculate tpr
res %>%
	map_dfr(
		~
			# Integrate
			.x$result$proportions %>%
			filter(level == 3) %>% 
			select(-draws, -rng, - .variable) %>% 
			unnest(proportions)%>%
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


# get reference from ARMET
ref =
	mix_base %>%
	inner_join(
		res[[1]]$result$internals$reference_filtered %>%
			filter(!`house keeping` ) %>%
			distinct(symbol)
	) %>%
	distinct(`Cell type category`, symbol,  `count normalised bayes`) %>%
	spread(`Cell type category`, `count normalised bayes`) %>%
	ttBulk::as_matrix(rownames = "symbol")

noiseles_test_ttBulk = function(mix, ref, method = "cibersort") {

	rr =
		mix %>%
		mutate(run = as.character(run)) %>%
		ttBulk::deconvolve_cellularity(run, symbol, `count mix`, reference = ref, method = method) %>%
		dplyr::select(run, days, alive, contains(method)) %>%
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
				~ .x %>%
					mutate(happened = as.integer(!alive)) %>%
					survival::coxph(survival::Surv(days, happened) ~ .value, .) %>%
					broom::tidy() %>%
					rename(alpha_2 = estimate) %>%
					mutate(		.lower_alpha2 = conf.low,	.upper_alpha2 = conf.high	) %>%
					dplyr::select(alpha_2, .lower_alpha2, .upper_alpha2, p.value)
			)
		) %>%
		unnest(cols = c(data, CI))

	list(mix = mix, result = rr)
}


# Same run with CIBERSORT

res_cibersort =
	mixes %>%
	map( ~ .x %>% noiseles_test_ttBulk(ref))

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
			#rename(run = sample) %>%
			#select(run, `Cell type category`, .value) %>%
			select(-alpha_2) %>%
			left_join(
				.x$mix %>% attr("proportions") %>% mutate(run = run %>% as.character)
			) %>%

			# Calculate
			mutate(inside = p > .lower_alpha2 & p < .upper_alpha2) %>%
			count(`Cell type category`, inside) %>%
			spread(inside, n) %>%
			replace(., is.na(.), 0) %>%
			mutate(tpr = `TRUE` / (`TRUE` + `FALSE`))
	) %>%
	ggplot(aes(tpr, x = `Cell type category`)) +
	geom_boxplot(outlier.shape = NA) +
	geom_jitter() +
	my_theme


# Same run with LLSR
res_llsr =
	which_is_up_down %>%
	map( ~ .x %>% noiseles_test_ttBulk(ref, method="llsr"))

res_llsr[[1]]$mix %>% attr("proportions") %>%
	mutate(run = run %>% as.character) %>%
	dplyr::select(-contains("alpha")) %>%
	left_join(res_llsr[[1]]$result) %>%
	ggplot(aes(x = p, y = .value, color = `Cell type category`)) +
	geom_abline(intercept = 0 , slope = 1) +
	geom_smooth(method = "lm") +
	geom_point()


# Create boxplot of errors
library(broom)
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
	mutate(error_logit_space = abs( gtools::logit(p) -  gtools::logit(.value))) %>%

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
	ggplot(aes(x = model, y = error_logit_space)) +
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


#------------------------------------#
# REGRESSION
#------------------------------------#

res_regression = readRDS("dev/test_student_noisless_regression.rds")

CI_to_ARMET = function(.data, CI){
	.data %>%
		map_dfr(
			~
				# Integrate
				.x$result %>%

				test_differential_composition(credible_interval = CI) %>%
				filter(level == 3) %>%
				left_join(
					.x$mix %>% attr("proportions") %>% dplyr::distinct(`Cell type category`, alpha_2),
					by = "Cell type category"
				) %>%
				drop_na  %>%

				# Calculate
				mutate(fp = alpha_2 == 0 & significant) %>%
				mutate(fn = alpha_2 != 0 & !significant)

		)
}

CI_to_others = function(.data, pvalue){
	.data %>%
		map_dfr(
			~
				# Integrate
				.x$result %>%
				dplyr::select(`Cell type category`, contains("alpha2"), 	`p.value`) %>%
				distinct() %>%
				left_join(
					.x$mix %>% attr("proportions") %>% distinct(`Cell type category`, alpha_2),
					by = "Cell type category"
				) %>%
				drop_na  %>%

				# Calculate
				mutate(fp = alpha_2 == 0 &
							 	`p.value` < pvalue) %>%
				mutate(fn = alpha_2 != 0 & `p.value` > pvalue)

		)
}

# Calculate fpr for regression
regression_ARMET =
	tibble(CI = c( seq(0.0, 0.85, 0.005), seq(0.85, 0.9999, 0.001))) %>%
	mutate(roc = future_map(
		CI,
		~ {
			x = CI_to_ARMET(res_regression, .x)

			x %>% filter(alpha_2 != 0) %>% group_by(alpha_2) %>%
				summarise(fnr = sum(fn) / n()) %>%
				mutate(
					fpr =
						x %>% filter(alpha_2 == 0) %>% select(-alpha_2) %>%
						summarise(fpr = sum(fp) / n()) %>% pull(1)
				)
		}
	))


regression_cibersort =
	tibble(CI = c( seq(0.0, 0.03, 0.001), seq(0.03, 0.9999, 0.001))) %>%
	mutate(roc = future_map(
		CI,
		~ {
			x = CI_to_others(res_cibersort, .x)

			x %>% filter(alpha_2 != 0) %>% group_by(alpha_2) %>%
				summarise(fnr = sum(fn) / n()) %>%
				mutate(
					fpr =
						x %>% filter(alpha_2 == 0) %>% select(-alpha_2) %>%
						summarise(fpr = sum(fp) / n()) %>% pull(1)
				)
		}
	))

regression_llsr =
	tibble(CI = c( seq(0.0, 0.03, 0.001), seq(0.03, 0.9999, 0.001))) %>%
	mutate(roc = future_map(
		CI,
		~ {
			x = CI_to_others(res_llsr, .x)

			x %>% filter(alpha_2 != 0) %>% group_by(alpha_2) %>%
				summarise(fnr = sum(fn) / n()) %>%
				mutate(
					fpr =
						x %>% filter(alpha_2 == 0) %>% select(-alpha_2) %>%
						summarise(fpr = sum(fp) / n()) %>% pull(1)
				)
		}
	))

regression_ARMET %>%
	mutate(method = "ARMET") %>%
	bind_rows(
		regression_cibersort %>%
			mutate(method = "cibersort")
	) %>%
	bind_rows(
		regression_llsr %>%
			mutate(method = "llsr")
	) %>%
	unnest(roc) %>%
	#bind_rows(tibble(CI = 1, alpha_2 = c(-4, 0.5, 1, 2), fnr = 1, fpr = 0)) %>%
	arrange(fnr) %>%
	mutate(tpr = 1-fnr) %>%
	ggplot(aes(x=fpr, y=tpr, color=method)) +
	geom_line() +
	scale_color_brewer(palette = "Set1") +
	facet_wrap(~ alpha_2, nrow = 1) +
	xlim(c(0,0.1)) +
	my_theme




CI_to_ARMET(res_regression, 0.95) %>%
	group_by(`Cell type category`, alpha_2) %>%
	summarise(fpr = sum(fp) / n(), fnr = sum(fn) / n()) %>%
	gather(which, rate, c("fpr", "fnr")) %>%
	ggplot(aes(y = rate, x = `Cell type category`, color = factor(alpha_2))) +
	geom_jitter() +
	facet_wrap( ~ which) +
	my_theme



# Calculate fpr for regression

regression_cibersort %>%
	group_by(`Cell type category`, alpha_2) %>%
	summarise(fpr = sum(fp) / n(), fnr = sum(fn) / n()) %>%
	gather(which, rate, c("fpr", "fnr"))%>%
	ggplot(aes(y = rate, x = `Cell type category`, color = factor( alpha_2))) +
	geom_jitter() +
	facet_wrap( ~ which) +
	my_theme


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


# Regression accuracy
CI_to_ARMET(res_regression, 0.95) %>% mutate(algorithm="ARMET") %>%
	bind_rows(	regression_cibersort %>% mutate(algorithm="cibersort")) %>%
	bind_rows(	regression_llsr %>% mutate(algorithm="llsr")) %>%
	group_by( alpha_2, algorithm) %>%
	summarise(fpr = sum(fp) / n(), fnr = sum(fn) / n()) %>%
	gather(which, rate, c("fpr", "fnr")) %>%
	ggplot(aes(y = rate, x = algorithm, color = factor(alpha_2))) +
	geom_point() +
	facet_grid( ~ which) +
	my_theme
