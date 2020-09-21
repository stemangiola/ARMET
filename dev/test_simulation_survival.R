# Create mix_base within NA - imputed values

library(tidyverse)
library(magrittr)
library(purrr)
library(furrr)
library(data.tree)
library(foreach)
library(ARMET)
plan(multiprocess, workers=28)
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


my_dir = "dev/test_simulation"


process_third_party = function(CI, file, input_file, slope){
	readRDS(file) %>% 
		
		# Add truth
		dplyr::select(.cell_type, estimate, 	p.value) %>%
		left_join(
			readRDS(input_file) %>% 
				attr("proportions") %>% 
				distinct(`Cell type category`, alpha_2),
			by = c(".cell_type" = "Cell type category"  )
		) %>%
		drop_na  %>%
		
		# Calculate
		mutate(fp = alpha_2 == 0 &	p.value < CI) %>%
		mutate(tp = alpha_2 != 0 &	p.value < CI)  %>%
		
		# Filter out accidental fp because of simplex
		filter(!(fp & (slope * estimate )<0)) 
}

process_ARMET = function(CI, file, input_file, slope){
	readRDS(file) %$%
		proportions %>% 
		#get_CI(CI) %>%
		
		mutate(.value_2 = map_dbl(
			draws_cens,
			~ .x %>% 
				filter(A == 2) %>%
				pull(.value ) %>%
				mean
		)) %>%
		
		mutate(prob_non_0_2 = map_dbl(
			draws_cens,
			~ .x %>% 
				filter(A == 2) %>% 
				ARMET:::draws_to_prob_non_zero()
		)) %>%
		
		# Add truth
		dplyr::select(`Cell type category`, .value_2, prob_non_0_2) %>%
		left_join(
			readRDS(input_file) %>% 
				attr("proportions") %>% 
				distinct(`Cell type category`, alpha_2),
			by = c("Cell type category"  )
		) %>%
		drop_na  %>%
		
		# Calculate
		mutate(fp = (alpha_2 == 0) &	abs(prob_non_0_2) > CI) %>%
		mutate(tp = (alpha_2 != 0) &	abs(prob_non_0_2) > CI) %>%
		
		# Filter out accidental fp because of simplex
		filter(!(fp & (slope * .value_2 )<0))  
	
}

roc_df = 
	
	dir(my_dir, full.names = TRUE) %>%
	enframe(value = "file", name = NULL) %>% 
	
	# format
	separate(file, as.character(1:7), sep = "__", remove = F) %>% 
	gather(which, value, -file) %>% 
	mutate(value = gsub("dev/test_simulation/", "", value)) %>% 
	mutate(value = gsub(".rds", "", value, fixed = T))  %>% 
	tidyr::extract(value, c("prefix", "value"), "([a-zA-z]+_)?(.+)") %>% 
	mutate(prefix = if_else(prefix == "", "in_out", gsub("_", "", prefix))) %>%
	select(-which) %>%
	spread(prefix, value) %>%
	nest(data = -c(foreignProp, run, S, slope, whichChanging)) %>%
	mutate(data = map(
		data,
		~.x %>% 
			filter(in_out=="output") %>% 
			mutate(input_file = 
						 	.x %>% 
						 	filter(in_out=="input") %>% 
						 	pull(file)
						)
	)) %>%
	unnest(data) %>%
	select(-`<NA>`) %>%
	mutate(slope = as.numeric(slope)) %>%
	
	# Add CI interval
	mutate(
		CI = 
			if_else(
				method=="ARMET",
				list(c( seq(0.0, 0.85, 0.05), seq(0.85, 0.9999, 0.001))),
				list(rev(1-c( seq(0.0, 0.85, 0.05), seq(0.85, 0.9999, 0.001),seq(0.9999, 0.99999999, 0.000001) )))
			)
	) %>%
	unnest(CI) %>%
	
	# Calculate
	mutate(fp_tp = future_pmap(
		list(CI, file, input_file, slope, method),
		~ ..5 %>% purrr::when(
			(.) == "ARMET" ~ process_ARMET(..1, ..2, ..3, ..4),
			~ process_third_party(..1, ..2, ..3, ..4)

	))) %>%
	unnest(fp_tp) %>%
	
	# Calculate ROC
	mutate(abs_slope = abs(slope)) %>%
	nest(data = -c(foreignProp, S, abs_slope, method, CI)) %>%
	mutate(real_negative = map_dbl(data, ~ .x %>% filter(alpha_2==0) %>% nrow)) %>%
	mutate(FP = map_dbl(data, ~ .x %>% filter(fp) %>% nrow)) %>%
	mutate(real_positive = map_dbl(data, ~ .x %>% filter(alpha_2!=0) %>% nrow )) %>%
	mutate(TP = map_dbl(data, ~ .x %>% filter(tp) %>% nrow) ) %>%
	mutate(TP_rate = TP/real_positive, FP_rate = FP/real_negative) %>%
	
	select(-data)

saveRDS(merged, "dev/test_simulation/roc_df.rds", compress = "xz")


# 
# roc_df %>%
# 	arrange(FP_rate, TP_rate) %>%
# 	ggplot(aes(x=FP_rate, y=TP_rate, color=method)) +
# 	geom_abline(intercept = 0, slope = 1, linetype="dotted", color="grey") +
# 	geom_line() +
# 	scale_color_brewer(palette = "Set1") +
# 	facet_wrap(~ abs_slope + foreignProp + S, nrow = 1) +
# 	coord_cartesian(xlim=c(0,0.08), ylim=c(0,1))
# 
# # merge results
# merged = 
# 	regression_ARMET %>%
# 	mutate(method="armet") %>%
# 	bind_rows(
# 		regression_cibersort %>%
# 			unnest(roc) %>%
# 			mutate(method = "cibersort") 
# 	) %>%
# 	arrange(FP_rate, TP_rate)
# 
# saveRDS(merged, "dev/test_simulation/merged.rda", compress = "xz")
# 
# (
# 	merged %>%
# 	ggplot(aes(x=FP_rate, y=TP_rate, color=method)) +
# 	geom_abline(intercept = 0, slope = 1, linetype="dotted", color="grey") +
# 	geom_line() +
# 	scale_color_brewer(palette = "Set1") +
# 	facet_wrap(~ abs_slope_run, nrow = 1) +
# 	coord_cartesian(xlim=c(0,0.08), ylim=c(0,1)) +
# 	my_theme
# )




