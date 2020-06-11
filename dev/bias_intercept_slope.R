# mean slope bias
# Is this due to prior on alpha?

library(tidyverse)
library(ARMET)
library(furrr)
library(tidybulk)
library(doParallel)
registerDoParallel(25)
# #foreach(i = dir("dev/armet_TCGA_may19_before_log_days/", pattern = "^armet_", full.names = T), .combine = bind_rows) %dopar% {
df = 
	foreach(i = dir("dev/", pattern = "^armet_", full.names = T) %>% grep("rda$", ., value = T), .combine = bind_rows) %do% {

	load(i)
	res %>% test_differential_composition() %>% select(.variable, community, level ,   `Cell type category`, C, .value_alpha1, .value_alpha2, zero) %>%
		mutate(file=i)
} 

df %>%
	#filter(level ==3) %>%
	#filter(`Cell type category` %in% c("t_CD8", "t_CD4", "macrophage") %>% `!`) %>%
	ggplot(aes( .value_alpha1, .value_alpha2-zero, label=`Cell type category`, color = factor(level))) +
	geom_point(aes()) +
	geom_smooth( method = "lm") +
	geom_abline(intercept = 0.042, slope = 0.8182)

# Regression full least square - PCA
#{ odregress(.$.value_alpha1, .$.value_alpha2-.$zero) } %$% coeff

# Test with simulated data and with real alpha intercept and 0 slope
load("/stornext/Home/data/allstaff/m/mangiola.s/PhD/deconvolution/ARMET/dev/armet_ACC.tcga.harmonized.counts.allgenes.rds.rda")

mix_base = ARMET:::get_noiseless_harmonised()
cell_types =  mix_base %>% filter(level ==3) %>% pull(`Cell type category`) %>% unique
alpha =  res %>% test_differential_composition() %>% filter(`Cell type category` %in% cell_types) %>% arrange(`Cell type category`) %>% select(.value_alpha1) %>% mutate(slope=0) %>% nanny::as_matrix()
X_df = res$proportions %>% select(proportions) %>% unnest(proportions) %>% select(sample, real_days = PFI.time.2, days = PFI.time.2, alive) %>% distinct %>% mutate(intercept = 1)

mix = mix_base %>% ARMET:::generate_mixture(X_df, alpha)

rr_no_buffer_beta1 =
	mix %>%
	mutate(`count mix` = as.integer(`count mix`), run = as.character(run)) %>%
	select(-level) %>%
	ARMET_tc(
		~ censored(days, alive),
		run, symbol, `count mix`,
		prior_survival_time = X_df$real_days %>% as.numeric
		#,
		#model = rstan::stan_model("~/PhD/deconvolution/ARMET/inst/stan/ARMET_tc_fix_hierarchical.stan", auto_write = F)
	) 

rr_no_buffer_beta2 = rr_no_buffer_beta1 %>%
	ARMET_tc_continue(2,model = rstan::stan_model("~/PhD/deconvolution/ARMET/inst/stan/ARMET_tc_fix_hierarchical.stan")) 

rr_no_buffer_beta3 = rr_no_buffer_beta2 %>%
	ARMET_tc_continue(3)
										#,model = rstan::stan_model("~/PhD/deconvolution/ARMET/inst/stan/ARMET_tc_fix_hierarchical.stan")) 

rr_no_buffer_beta4 = rr_no_buffer_beta3 %>%
	ARMET_tc_continue(4)
										#,model = rstan::stan_model("~/PhD/deconvolution/ARMET/inst/stan/ARMET_tc_fix_hierarchical.stan")) 


rr_no_buffer %>%
	test_differential_composition()  %>%
	filter(level ==3) %>% 
	filter(`Cell type category` %in% c("t_CD8", "t_CD4", "macrophage") %>% `!`) %>%
	ggplot(aes( .value_alpha1, .value_alpha2-zero, label=`Cell type category`)) +
	geom_point(aes(color = factor(level))) +
	geom_smooth( method = "lm") + 
	geom_abline(intercept = 0.042, slope = 0.8182)

(rr_no_buffer_beta2$proportions %>%
		unnest(draws) %>%
		filter(A == 2) %>%
		ggplot(aes(.value, color=`Cell type category`)) +
		geom_density() +
		facet_wrap(~.variable, scale="free_y")
) %>% plotly::ggplotly()

rr %>% plot_scatter() + scale_x_log10() 
