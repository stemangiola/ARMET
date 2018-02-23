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
	
	# Get ref of the current level
	ref = ref %>%
		dplyr::mutate(ct = as.character(ct)) %>%
		dplyr::left_join(get_map_foreground_background(my_tree, ct), by="ct") %>%
		dplyr::rename(ct_ = ct) %>%
		dplyr::filter(gene %in% 
										get_node_label_level_specfic(
											node_from_name(my_tree,  ct), 
											label = "markers",
											start_level = 1, 
											stop_level = 1
										)
		) %>%
		dplyr::select(-ct_) %>%
		dplyr::rename(ct = ancestor) %>%
		dplyr::mutate_if(is.character, as.factor)
	
	
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
		all(
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
				dplyr::group_by(gene, ct) %>%
				dplyr::summarise(value = median(value)) %>%
				dplyr::ungroup() %>%
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
			dplyr::select(-sample),
		
		x_genes = fg %>% 
			dplyr::pull(gene) %>%
			levels() %>%
			tibble::as_tibble() %>%
			dplyr::rename(gene = value),
		
		x = fg %>% 
			dplyr::select(gene, ct, value) %>%
			dplyr::group_by(gene, ct) %>%
			dplyr::summarise(value = median(value)) %>%
			dplyr::ungroup() %>%
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
		is_mix_microarray =  as.numeric(is_mix_microarray)
		
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

	# Add hypothesis testing
	node = 
		if(!is.null(cov_to_test))
			add_data_to_tree_from_table(
				node,
				dirReg_test(
					fit, 
					my_design, 
					cov_to_test, 
					names_groups = levels(fg$ct) 
				)$stats,
				"stats",
				append = F
			)
		else node

	
	# Save output
	if(save_report) save(fit, file=sprintf("%s/%s_fit.RData", output_dir, ct))
	# p = rstan::traceplot(fit, pars=c("alpha"), inc_warmup=F)
	# 	
	#if(save_report) ggplot2::ggsave(sprintf("%s_chains.png", ct), p)
	if(save_report)	write.csv(proportions, sprintf("%s/%s_composition.csv", output_dir, ct))
	
	
	list(proportions = proportions, node = node)
}
