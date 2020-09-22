# Create mix_base within NA - imputed values

# # Build HPC list	
# tibble(block = c("CORES=12\nMEMORY=30024")) %>%
# 	bind_rows(
# 		expand.grid(c(0.5, 1, 2, -0.5, -1, -2),1:16) %>%
# 			mutate(file = sprintf("dev/test_simulation/test_noisless_survival_regression_%s_%s.rds:", Var1, Var2)) %>%
# 			mutate(command = sprintf("\n\tsh dev/test_noise_survival_standalone.sh %s %s", Var1, Var2)) %>%
# 			unite("block", c(file, command), sep="") %>%
# 			select(block)
# 	) %>%
# 	pull(block) %>%
# 	writeLines("dev/test_noise_survival_standalone.makeflow")

#~/third_party_sofware/cctools-7.1.5-x86_64-centos6/bin/makeflow  -T slurm -B '--time=10:60:00' --do-not-save-failed-output dev/test_noise_survival_standalone.makeflow
# ~/unix3XX/third_party_sofware/cctools-7.1.5-x86_64-centos6/bin/makeflow  -T torque -B '-l walltime=10:60:00' --do-not-save-failed-output dev/test_noise_survival_standalone.makeflow

args = commandArgs(trailingOnly=TRUE)
slope = as.numeric(args[1])
run = as.integer(args[2])
foreign_prop = 0.5
n_samples = 60

library(tidyverse)
library(magrittr)
library(purrr)
library(furrr)
library(data.tree)
library(foreach)
library(ARMET)
library(nanny)
library(tidybulk)

# plan(multisession, workers=10)

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
		filter(PFI.2 == 1 & !is.na(PFI.time.2) & PFI.time.2 > 0) %>%
		select(real_days = PFI.time.2 ) %>%
		mutate(real_days = real_days %>% scale(center = F) %>% as.numeric) %>%
		sample_n(S) %>%
		mutate(sample = sprintf("S%s", 1:n())) %>%
		mutate(alive = sample(0:1, n(), replace = T)) %>%
		mutate(days = ifelse(alive==1, real_days/2, real_days) ) %>%
		mutate(intercept = 1)
}

generate_mixture = function(.data, X_df, alpha, foreign_prop = 0) {
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
	
	# Regress on the log days
	X = X_df %>% mutate(real_days = log(real_days)) %>% select(intercept, real_days) %>% nanny::as_matrix()
	
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
	
	ct_changing = alpha_df %>% filter(alpha_2 != 0) %>% pull(`Cell type category`)
	
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
	
	# Add fold increase decrease
	fold_change = 
		matrix(c(rep(1, 2), c(0, max(X,2))), ncol = 2)  %*% t(alpha) %>%
		apply(1, softmax) %>%
		t %>%
		`*` (40) %>%
		apply(1, softmax) %>%
		.[ct_names == ct_changing,] %>%
		{	max(.) / min(.)	} %>%
		{ slope = alpha[,2][ alpha[,2]!=0]; ifelse(slope<0, -(.), (.)) }
	
	# Add counts
	dirichlet_source =
		cell_type_proportions %>%
		left_join(.data, by = c("Cell type category", "sample"))
	
	# Add foreign sample
	neural_sample = 
		readRDS("dev/ENCODE.rds") %>% 
		filter(grepl("neur", `Cell type`)) %>% 
		group_by(symbol) %>% slice(1) %>% ungroup() %>%
		select(symbol, count_neuro = count)
	
	# Make mix
	dirichlet_source %>%
		mutate(c = `count scaled bayes` * p) %>%
		group_by(run, symbol) %>%
		summarise(`count mix` = c %>% sum) %>%
		ungroup %>%
		
		left_join(dirichlet_source %>% nanny::subset(run) ) %>%
		
		mutate(fold_change = fold_change) %>%
		
		# Add neuron
		left_join(neural_sample, by = "symbol") %>%
		mutate(prop_neural = foreign_prop) %>%
		mutate(`count mix` = (`count mix` * (1-prop_neural))+(count_neuro*prop_neural)) %>%
		
		# Add proportions
		add_attr(cell_type_proportions, "proportions") 
	

}

noise_test = function(mix_base, slope, which_changing, S, foreign_prop = 0) {
	
	cell_types =  mix_base %>% filter(level ==3) %>% pull(`Cell type category`) %>% unique
	
	alpha = get_alpha(slope, which_changing, cell_types)
	X_df = get_survival_X(S)
	
	mix = mix_base %>% generate_mixture(X_df, alpha, foreign_prop = foreign_prop)
	
	rr =
		mix %>%
		mutate(`count mix` = as.integer(`count mix`), run = as.character(run)) %>%
		select(-level) %>%
		ARMET_tc(
			~ censored(days, alive),
			run, symbol, `count mix`,
			prior_survival_time = X_df$real_days %>% as.numeric, iterations = 2000, sampling_iterations = 500
		)  %>%
		ARMET_tc_continue(2) %>%
		ARMET_tc_continue(3)
	
	rr$proportions = rr %>% add_cox_test("days", "alive")
	
	list(mix = mix, result = rr)
	#%>% saveRDS(sprintf("dev/test_student_noisless_%s", which_is_up_down %>% paste(collapse="_")))
}

get_noiseless_harmonised = function(){
	
	mix_base_unharmonized = readRDS("dev/mix_base_noiseless.RDS")
	
	my_markers =
		ARMET::ARMET_ref %>%
		
		left_join(ARMET::n_markers, by = c("ct1", "ct2")) %>%
		filter_reference(
			mix_base_unharmonized %>%
				filter(level == 3) %>%
				distinct(`Cell type category`, symbol, `count normalised bayes`) ,
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

noise_test_cibersort = function(mix_base, slope, which_changing, S, method="cibersort", foreign_prop = 0) {
	
	cell_types =  mix_base %>% filter(level ==3) %>% pull(`Cell type category`) %>% unique
	
	alpha = get_alpha(slope, which_changing, cell_types)
	X_df = get_survival_X(S)
	
	mix = mix_base %>% generate_mixture(X_df, alpha, foreign_prop = foreign_prop)
	
	rr =
		mix %>%
		mutate(`count mix` = as.integer(`count mix`), run = as.character(run)) %>%
		select(-level) %>%
		mutate(dead = !alive) %>%
		tidybulk::test_differential_cellularity(
			survival::Surv(days, dead) ~ .,
			run, symbol, `count mix`,
			reference = get_noiseless_harmonised() %>%
				distinct(`Cell type category`, symbol, `count normalised bayes`) %>%
				spread(`Cell type category`, `count normalised bayes`) %>%
				nanny::as_matrix(rownames = symbol) %>%
				as.data.frame,
			method = method
		)

	list(mix = mix, result = rr)

}

mix_base = readRDS("dev/mix_base.RDS") %>% filter(level==3)

noise_test(mix_base, slope, run, n_samples, foreign_prop = foreign_prop) %>%
	saveRDS(sprintf("dev/test_simulation/test_noise_survival_regression_%s_%s.rds", slope, run), compress = "gzip")

which_is_up_down = 1:16

# future_map(which_is_up_down,  ~ noise_test_cibersort(mix_base, 2, .x, n_samples, method="cibersort", foreign_prop = foreign_prop)) %>%
# 	c(future_map(which_is_up_down,  ~ noise_test_cibersort(mix_base, -2, .x, n_samples, method="cibersort", foreign_prop = foreign_prop))) %>%
# 	saveRDS("dev/test_simulation/test_noise_survival_regression_2_slope_cibersort.rds")
# 
# 
# future_map(which_is_up_down,  ~ noise_test_cibersort(mix_base, 1, .x, n_samples, method="cibersort", foreign_prop = foreign_prop)) %>%
# 	c(future_map(which_is_up_down,  ~ noise_test_cibersort(mix_base, -1, .x, n_samples, method="cibersort", foreign_prop = foreign_prop))) %>%
# 	saveRDS("dev/test_simulation/test_noise_survival_regression_1_slope_cibersort.rds")
# 
# 
# future_map(which_is_up_down,  ~ noise_test_cibersort(mix_base, 0.5, .x, n_samples, method="cibersort", foreign_prop = foreign_prop)) %>%
# 	c(future_map(which_is_up_down,  ~ noise_test_cibersort(mix_base, -0.5, .x, n_samples, method="cibersort", foreign_prop = foreign_prop))) %>%
# 	saveRDS("dev/test_simulation/test_noise_survival_regression_0.5_slope_cibersort.rds")
# 



# 
# future_map(which_is_up_down,  ~ noise_test_cibersort(mix_base, 2, .x, n_samples, method="llsr", foreign_prop = foreign_prop)) %>%
# 	c(future_map(which_is_up_down,  ~ noise_test_cibersort(mix_base, -2, .x, n_samples, method="llsr", foreign_prop = foreign_prop))) %>%
# 	saveRDS("dev/test_simulation/test_noise_survival_regression_2_slope_llsr.rds")
# 
# 
# future_map(which_is_up_down,  ~ noise_test_cibersort(mix_base, 1, .x, n_samples, method="llsr", foreign_prop = foreign_prop)) %>%
# 	c(future_map(which_is_up_down,  ~ noise_test_cibersort(mix_base, -1, .x, n_samples, method="llsr", foreign_prop = foreign_prop))) %>%
# 	saveRDS("dev/test_simulation/test_noise_survival_regression_1_slope_llsr.rds")
# 
# 
# future_map(which_is_up_down,  ~ noise_test_cibersort(mix_base, 0.5, .x, n_samples, method="llsr", foreign_prop = foreign_prop)) %>%
# 	c(future_map(which_is_up_down,  ~ noise_test_cibersort(mix_base, -0.5, .x, n_samples, method="llsr", foreign_prop = foreign_prop))) %>%
# 	saveRDS("dev/test_simulation/test_noise_survival_regression_0.5_slope_llsr.rds")
# 

# 
# 
# 
# future_map(which_is_up_down,  ~ noise_test_cibersort(mix_base, 2, .x, n_samples, method="epic", foreign_prop = foreign_prop)) %>%
# 	c(future_map(which_is_up_down,  ~ noise_test_cibersort(mix_base, -2, .x, n_samples, method="epic", foreign_prop = foreign_prop))) %>%
# 	saveRDS("dev/test_simulation/test_noise_survival_regression_2_slope_epic.rds")
# 
# 
# future_map(which_is_up_down,  ~ noise_test_cibersort(mix_base, 1, .x, n_samples, method="epic", foreign_prop = foreign_prop)) %>%
# 	c(future_map(which_is_up_down,  ~ noise_test_cibersort(mix_base, -1, .x, n_samples, method="epic", foreign_prop = foreign_prop))) %>%
# 	saveRDS("dev/test_simulation/test_noise_survival_regression_1_slope_epic.rds")
# 
# 
# future_map(which_is_up_down,  ~ noise_test_cibersort(mix_base, 0.5, .x, n_samples, method="epic", foreign_prop = foreign_prop)) %>%
# 	c(future_map(which_is_up_down,  ~ noise_test_cibersort(mix_base, -0.5, .x, n_samples, method="epic", foreign_prop = foreign_prop))) %>%
# 	saveRDS("dev/test_simulation/test_noise_survival_regression_0.5_slope_epic.rds")