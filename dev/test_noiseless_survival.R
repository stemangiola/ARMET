# Create mix_base withn NA - imputed values

library(tidyverse)
library(magrittr)
library(purrr)
library(furrr)
library(data.tree)
library(foreach)
library(ARMET)
# plan(multiprocess)
library(doParallel)
registerDoParallel(cores=20)
library(iterators)
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

#
# which_is_up_down %>%
# 	map( ~ .x %>% noiseles_test(do_regression = T)) %>%
# 	saveRDS("dev/test_student_noisless_regression.rds")

#------------------------------------#
# REGRESSION
#------------------------------------#

dir("dev", "test_noisless_survival_regression", full.names = T) %>%
	map_dfr(~ .x %>% readRDS %>% mutate(file=.x))

res_regression = readRDS("dev/test_noisless_survival_regression.rds")

CI_to_ARMET = function(.data, CI){
	
	foreach(x = .data, i=icount(), .combine = bind_rows) %dopar% {
		x$result %>%
			
			test_differential_composition(credible_interval = CI) %>%
			#filter(level == 3) %>%
			inner_join(
				x$mix %>% attr("proportions") %>% dplyr::distinct(`Cell type category`, alpha_2) %>%
					
					# Adapt to the detection of cell imbalance rather than absolute change
					{
						my_ct =(.) %>% filter(alpha_2 != 0) %>% pull(`Cell type category`)
						my_ancestors = ARMET::tree %>% ARMET:::ToDataFrameTypeColFull() %>% filter(level_4 == my_ct) %>% as.character
						my_cousin = ARMET::tree %>% ARMET:::ToDataFrameTypeColFull() %>% filter(level_3 == my_ancestors[3]) %>% pull(level_4) 
						my_all = purrr::when(length(my_cousin) > 2 ~ my_ancestors, ~ c(my_ancestors, my_cousin)) %>% unique
						my_slope = (.) %>% filter(alpha_2!=0) %>% pull(alpha_2)
						(.) %>% mutate(alpha_2 = if_else(`Cell type category`%in% my_all, my_slope, alpha_2))
					}	,
				by = "Cell type category"
			) %>%
		#	drop_na  %>%
			
			# Calculate
			mutate(fp = alpha_2 == 0 & significant) %>%
			mutate(tp = alpha_2 != 0 & significant) %>%
			
			mutate(run = i)
		
	}
	

}

# Calculate fpr for regression
regression_ARMET =
	tibble(CI = c( seq(0.0, 0.85, 0.05), seq(0.85, 0.9999, 0.001))) %>%
	#slice(1, 10, 30, 100) %>%
	mutate(roc = map(
		CI,
		~ {
			cat(".")
			x = CI_to_ARMET(res_regression, .x)

			slope_run = x %>% distinct(alpha_2, run) %>% filter(alpha_2 != 0) %>% rename(slope_run = alpha_2)
			
			x %>%
				left_join(slope_run, by="run") %>%
				nest(data = -slope_run) %>%
				mutate(real_negative = map_dbl(data, ~ .x %>% filter(alpha_2==0) %>% nrow)) %>%
				mutate(FP = map_dbl(data, ~ .x %>% filter(fp) %>% nrow)) %>%
				mutate(real_positive = map_dbl(data, ~ .x %>% filter(alpha_2!=0) %>% nrow )) %>%
				mutate(TP = map_dbl(data, ~ .x %>% filter(tp) %>% nrow) ) %>%
				mutate(TP_rate = TP/real_positive, FP_rate = FP/real_negative) %>%
			
				select(-data)
				
		}
	))


(
	regression_ARMET %>%
	mutate(method = "ARMET") %>%
	unnest(roc) %>%
	arrange(FP_rate) %>%
	ggplot(aes(x=FP_rate, y=TP_rate)) +
	geom_abline(intercept = 0, slope = 1, linetype="dotted", color="grey") +
	geom_line() +
	scale_color_brewer(palette = "Set1") +
	facet_wrap(~ slope_run, nrow = 1) +
	xlim(c(0,1)) +
		ylim(c(0,1)) +
	my_theme
) %>%
	ggsave(
		"dev/test_noisless_survival_ROC.pdf",
		plot = .,
		useDingbats=FALSE,
		units = c("mm"),
		width = 183 ,
		limitsize = FALSE
	)
