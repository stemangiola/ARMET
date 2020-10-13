library(tidyverse)
library(ARMET)
library(magrittr)
library(tidybulk)
library(tidyseurat)

args = commandArgs(trailingOnly=TRUE)
slope = as.numeric(args[1])
foreign_prop = as.numeric(args[2])
S = as.integer(args[3])
which_changing = as.integer(args[4])
run = as.integer(args[5])
output_file = args[6]

set.seed(34543*run)

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
		dplyr::sample_n(S) %>%
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
	
	# Avoid map error
	.my_data = .data
	
	# Regress on the log days
	X = X_df %>% mutate(real_days = log(real_days)) %>% select(intercept, real_days) %>% nanny::as_matrix()
	
	sample_list = .data %>% tidyseurat::pull(sample) %>% unique %>% sample(size = nrow(X), replace = TRUE)
	

	
	samples_run_df = 
		.data %>%
		tidyseurat::count(sample, name = "cell_total") %>%
		dplyr::sample_n(30, replace = T) %>%
		tibble::rowid_to_column(var = "run")
	
	ct_names = .data %>% tidyseurat::distinct(cell_type_curated) %>% pull(1)
	
	alpha_df = alpha %>% as.data.frame %>% setNames(sprintf("alpha_%s", 1:2)) %>% mutate(cell_type_curated  = ct_names)
	
	ct_changing = alpha_df %>% filter(alpha_2 != 0) %>% pull(cell_type_curated)
	
	cell_type_proportions =
		# Choose samples
		samples_run_df %>%
		
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
				gather(cell_type_curated, alpha, -run)
		) %>%
		
		# Add X
		left_join(X_df %>% select(-sample) %>% mutate(run = 1:n())) %>%
		
		# Add alpha
		left_join(alpha_df) %>%
		
		group_by(run) %>%
		mutate(p = gtools::rdirichlet(1, alpha) %>% as.numeric()) %>%
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
		mutate(cell_count = as.integer(cell_total*p)) %>%
		mutate(transcriptome = pmap(
			list(sample, cell_type_curated, cell_count),
			~ .my_data %>% 
				
				# Subset sample/cell-type
				tidyseurat::filter(sample == ..1 & cell_type_curated == ..2) %>%
				
				# take n cells based on proportion
				tidyseurat::sample_n(..3, replace = TRUE) %>%
				
				# Get abundance sum
				`@` (assays) %$% 
				RNA %>%
				`@` (counts) %>%
				as.matrix() %>%
				rowSums %>% 
				enframe(name = "transcript", value = "count")
		))

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

mix_base = readRDS("~/PhD/deconvolution/ARMET/dev/test_simulation_singleCell_makeflow_pipeline/PBMC_curate_subset.rds") 

cell_types =  mix_base %>% pull(cell_type_curated) %>% unique
alpha = get_alpha(slope, which_changing, cell_types)
X_df = get_survival_X(S)

mix_base %>%
	generate_mixture(X_df, alpha, foreign_prop = foreign_prop) %>%
	
	# Parse
	mutate(`count mix` = as.integer(`count mix`), run = as.character(run)) %>%
	select(-level) %>%
	rename(sample = run) %>%
	
	# Save
	saveRDS(output_file, compress = "xz")
	
