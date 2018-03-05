# ARMET_tc_coreAlg.R

ARMET_tc_coreAlg = function(
	obj.in, 
	node, 
	real_prop_obj=NULL, 
	is_test=F){
	
	# Get ct
	ct = node$name
	
	# Parse input
	mix =                   obj.in$mix
	ref =                   obj.in$ref
	my_design=              obj.in$my_design
	cov_to_test =           obj.in$cov_to_test
	fully_bayesian =        obj.in$fully_bayesian
	observed_prop =         obj.in$observed_prop
	ct_to_omit =            obj.in$ct_to_omit
	my_tree =               obj.in$my_tree
	is_mix_microarray =     obj.in$is_mix_microarray	
	output_dir =            obj.in$output_dir
	sigma_hyper_sd =        obj.in$sigma_hyper_sd
	phi_hyper_sd =          obj.in$phi_hyper_sd
	alpha_hyper_value =     obj.in$alpha_hyper_value
	save_report =           obj.in$save_report
	omit_regression =       obj.in$omit_regression
	do_debug =              obj.in$do_debug

	# Get ref of the current level
	ref = 
		ref %>%
		# Get cell type just for this recursive peers and descendants
		dplyr::mutate(ct = as.character(ct)) %>%
		dplyr::left_join(
			get_map_foreground_background(tree, ct) %>% 
				dplyr::select(-ct) %>%
				dplyr::rename(ct=ancestor) %>%
				dplyr::distinct(), 
			by="ct"
		) %>%
		# Filter non used cell types
		dplyr::filter(!is.na(variable)) %>%
		# Filter non used genes
		dplyr::filter(
			gene %in% 
			get_node_label_level_specfic(
				get_node_from_name(my_tree,  !!ct), 
				label = "markers",
				start_level = 1, 
				stop_level = 1
			)
		) %>%
		dplyr::mutate_if(is.character, as.factor) %>%
		droplevels() 
	
	# Print stats on ref
	get_stats_on_ref(ref,my_tree) %>% 
		dplyr::filter(ct %in% unique(ref$ct)) 
	
	# filter mix and add theta value
	mix = mix %>%
		dplyr::group_by(sample) %>%
		dplyr::mutate(
			theta=
				length(which(value>0)) / 
				n()
		) %>%
		dplyr::ungroup() %>%
		dplyr::filter(gene %in% unique(ref$gene))
	
	fg = ref %>% 
		dplyr::filter(variable=="main") %>%
		droplevels()
	
	# Setup background
	bg = 	ref %>% 
		dplyr::filter(variable=="background") %>%
		droplevels()
	
	# Garbage collection
	rm(ref)
	gc()
	
	# Get the probability table of the previous run
	ancestor_run_prop_table = get_last_existing_leaves_with_annotation(my_tree)
	
	# Sanity check
	if(
		any(
			ancestor_run_prop_table %>% 
			dplyr::group_by(sample) %>% 
			dplyr::summarise(tot=sum(absolute_proportion)) %>%
			dplyr::pull(tot) != 1 
		)
	) {
		writeLines("ARMET: The absolute proportions are supposed to sum to 1 for each sample")
		print(
			ancestor_run_prop_table %>% 
				dplyr::group_by(sample) %>% 
				dplyr::summarise(tot=sum(absolute_proportion))
		)
		#stop()
	}
	
	# Get the probability table of my cell type
	fg_prop = ancestor_run_prop_table %>%	
		dplyr::filter(ct==!!ct) %>%
		droplevels()
	
	# Get the probability table of the background
	bg_prop = ancestor_run_prop_table %>% 
		dplyr::filter(ct != !!ct) %>%
		droplevels()
	
	# Prepare input for full bayesian 
	# e.obj = prepare_input(fg, my_tree)
	# bg.obj = prepare_input(bg, my_tree)
	
	# Calculate the value of the genes for background 
	# !!rlang::sym(ct) -> for variable name without quotes

	y_hat_background = 
		as.matrix(
			bg_prop %>%
			{ 
				if(nrow(.)==0) 
					dplyr::bind_rows(
						fg_prop %>%
							dplyr::mutate(
								ct = factor("bg"), 
								relative_proportion = 0, 
								absolute_proportion = 0
							)
					)
				else . 
			} %>%
				dplyr::select(-relative_proportion) %>%
				tidyr::spread(ct, absolute_proportion) %>%
				dplyr::select(-sample)
		) %*% 
		as.matrix(
			bg %>% 
			{ 
				if(nrow(.)==0) 
					# Take fg as template
					fg %>%
					dplyr::distinct(gene) %>%
					dplyr::mutate(
						sample="s0",
						value = 0,
						ct = "bg",
						variable = "background"
					)
				else . 
			} %>%
				dplyr::mutate_if(is.character, as.factor) %>%
				dplyr::select(gene, ct, value) %>%
				tidyr::spread(gene, value) %>%
				dplyr::select(-ct) %>%
				dplyr::mutate_all(dplyr::funs(ifelse(is.na(.), 0, .)))
		) %>%
		tibble::as_tibble() %>%
		dplyr::mutate(sample = levels(mix$sample)) %>%
		dplyr::select(sample, dplyr::everything())

	# Create input object for the model
	model.in = list(
		
		G = fg %>% 
			dplyr::distinct(gene) %>% 
			nrow(),
		
		S = mix %>% 
			dplyr::distinct(sample) %>% 
			nrow(),
		
		P = fg %>% 
			dplyr::distinct(ct) %>% 
			nrow(),
		
		R = my_design %>% 
			dplyr::select(-sample) %>%
			ncol(),
		
		y = mix %>% 
			dplyr::select(gene, sample, value) %>% 
			tidyr::spread(gene, value) %>%
			dplyr::arrange(sample) %>%
			dplyr::select(-sample), 
		
		X = my_design %>% 
			dplyr::arrange(sample) %>%
			dplyr::select(-sample) %>%
			# transform factors into numeric safely
			dplyr::mutate_if(is.factor, as.character) %>%
			dplyr::mutate_if(is.character, as.numeric),
		
		x_genes = fg %>% 
			dplyr::pull(gene) %>%
			levels() %>%
			tibble::as_tibble() %>%
			dplyr::rename(gene = value),
		
		x = fg %>% 
			dplyr::select(gene, ct, value) %>%
			tidyr::spread(ct, value) %>%
			dplyr::mutate_if(is.numeric, dplyr::funs(ifelse(is.na(.), 0, .))) %>%
			dplyr::select(-gene),
		
		y_hat_background = y_hat_background %>%
			dplyr::select(-sample),
		
		p_target = fg_prop %>% 
			dplyr::pull(absolute_proportion) %>% 
			as.array(),
		
		theta = mix %>%
			dplyr::distinct(sample, theta) %>%
			dplyr::pull(theta) %>% 
			as.array(),
		
		sigma_hyper_sd =     sigma_hyper_sd,
		phi_hyper_sd =       phi_hyper_sd,
		alpha_hyper_value =  alpha_hyper_value,
		is_mix_microarray =  as.numeric(is_mix_microarray),
		omit_regression =    as.numeric(omit_regression)
		
	)
	
	# Save the raw model resul for debugging
	if(save_report) save(model.in, file=sprintf("%s/%s_model_in.RData", output_dir, ct))
	
	# Choose model
	model = switch(fully_bayesian + 1, stanmodels$ARMET_tcFix_recursive, stanmodels$ARMET_tc_recursive) 

	# Run model
	fit = 
		rstan::sampling(
			model,
			data=                             model.in,
			iter=                             1000 ,
			#control =                         list(adapt_delta = 0.99, stepsize = 0.01, max_treedepth =15),
			cores=4
		)
	
	#if(do_debug) browser()
	
	# Parse results
	proportions =  
		parse_summary_vector_in_2D(apply( as.matrix(fit, pars = "beta"), 2, mean)) %>%
		tibble::as_tibble() %>%
		setNames(levels(mix$sample)) %>%
		dplyr::mutate(ct = levels(fg$ct)) %>%
		tidyr::gather(sample, relative_proportion, -ct) %>%
		dplyr::mutate_if(is.character, as.factor)
	
	# Add info to the node
	node = add_data_to_tree_from_table(node, proportions, "relative_proportion", append = T)
	
	# Set up background trees
	node = add_absolute_proportions_to_tree(node)

	test = switch (
		is.null(cov_to_test) + 1,
		dirReg_test(
			fit, 
			my_design, 
			cov_to_test, 
			names_groups = levels(fg$ct) 
		),
		NULL
	) 
	
	# Add hypothesis testing
	node = switch(
		is.null(cov_to_test) + 1,
		add_data_to_tree_from_table(
			node,
			test$stats,
			"stats",
			append = F
		),
		node
	)

	# Add hypothesis testing
	node = switch(
		is.null(cov_to_test) + 1,
		add_info_to_tree(
			node, 
			ct, 
			"plot_props", 
			test$plot_props,
			append = F
		),
		node
	)
	
	# Add hypothesis testing
	node = switch(
		is.null(cov_to_test) + 1,
		add_info_to_tree(
			node, 
			ct, 
			"plot_coef_ang", 
			test$plot_coef_ang,
			append = F
		),
		node
	)
	
	node = switch(
		is.null(cov_to_test) + 1,
		add_info_to_tree(
			node, 
			ct, 
			"estimate_prop_with_uncertanties", 
			parse_fit_for_quantiles(fit, 0.5, "mean", my_design, levels(fg$ct)) %>%
				dplyr::left_join(	parse_fit_for_quantiles(fit, 0.95, "upper", my_design, levels(fg$ct)), by=c("sample", "ct")) %>%
				dplyr::left_join(	parse_fit_for_quantiles(fit, 0.05, "lower", my_design, levels(fg$ct)), by=c("sample", "ct")) %>%
				dplyr::left_join(	my_design %>% dplyr::select(-`(Intercept)`), by="sample"),
			append = F
		),
		node
	)
	
	# Add hypothesis testing
	node = switch(
		is.null(cov_to_test) + 1,
		add_info_to_tree(
			node, 
			ct, 
			"coef_ang_posterior_adj", 
			test$coef_ang_posterior_adj,
			append = F
		),
		node
	)
	
	# Save output
	if(save_report) save(fit, file=sprintf("%s/%s_fit.RData", output_dir, ct))
	# p = rstan::traceplot(fit, pars=c("alpha"), inc_warmup=F)
	# 	
	#if(save_report) ggplot2::ggsave(sprintf("%s_chains.png", ct), p)
	if(save_report)	write.csv(proportions, sprintf("%s/%s_composition.csv", output_dir, ct))
	
	
	list(proportions = proportions, node = node)
}
