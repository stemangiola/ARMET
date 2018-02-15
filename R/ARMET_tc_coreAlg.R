# ARMET_tc_coreAlg.R

ARMET_tc_coreAlg = function(
	obj.in, 
	ct, 
	real_prop_obj=NULL, 
	is_test=F){

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
	bg_tree =               obj.in$bg_tree
	sigma_hyper_sd =        obj.in$sigma_hyper_sd
	phi_hyper_sd =          obj.in$phi_hyper_sd
	alpha_hyper_value =     obj.in$alpha_hyper_value
	save_report =           obj.in$save_report
	
	# Get ref of the current level
	ref = ref %>%
		dplyr:::mutate(ct = as.character(ct)) %>%
		dplyr:::left_join(get_map_foreground_background(tree, ct), by="ct") %>%
		dplyr:::rename(ct_ = ct) %>%
		dplyr:::filter(gene %in% 
					 	get_node_label_level_specfic(
					 		node_from_name(tree,  ct), 
					 		label = "markers",
					 		start_level = 1, 
					 		stop_level = 1
					 	)
					) %>%
		dplyr:::select(-ct_) %>%
		dplyr:::rename(ct = ancestor) %>%
		dplyr:::mutate_if(is.character, as.factor)
		
	
	# Print stats on ref
	get_stats_on_ref(ref,tree) %>% 
		dplyr:::filter(ct %in% unique(ref$ct)) 
	
	# filter mix and add theta value
	mix = mix %>%
		dplyr:::group_by(sample) %>%
		dplyr:::mutate(
			theta=
				length(which(value>0)) / 
				n()
		) %>%
		ungroup() %>%
		filter(gene %in% unique(ref$gene))
	
	fg = ref %>% 
		dplyr:::filter(variable=="main") %>%
		droplevels()
	
	# Prepare input for full bayesian 
	e.obj = prepare_input(fg, tree)
	
	# Setup background
	bg = 	ref %>% 
		dplyr:::filter(variable=="background") %>%
		droplevels()
	
	# Set up background trees
	bg_trees = divide_trees_proportion_across_many_trees(bg_tree)
	bg_trees = lapply(bg_trees, add_absolute_proportions_to_tree)
	
	# Set up the proportions of the background
	beta_ancestor_run = as.data.frame(lapply(bg_trees, get_last_existing_leaves_with_annotation)) %>%
		tibble:::as_tibble() %>%
		{ 
			if(nrow(.)==0) 
				tibble:::tibble(rep(1, ncol(mix)), rep(0, ncol(mix))) %>%
					stats:::setNames(., c(ct, "bg"))  %>%
					dplyr:::mutate(
						sample=
							mix %>% 
								dplyr:::pull(sample) %>% 
								levels(.)
						)
			else . 
		} %>%
		dplyr:::select(sample, everything()) %>%
		tidyr:::gather(ct, value, -1) %>%
		dplyr:::mutate_if(is.character, as.factor)
	
	# Background dataframe for full bayesian
	bg.obj = prepare_input(bg, tree)
	
	# Calculate the value of the genes for background
	# !!rlang::sym(ct)
	y_hat_background = 
		as.matrix(
			beta_ancestor_run %>% 
				dplyr:::filter(ct != !!ct) %>%
				tidyr:::spread(ct, value) %>%
				dplyr:::select(-sample)
		) %*% 
		as.matrix(
			bg %>% 
				{ 
					if(nrow(.)==0) 
						ref %>%
						dplyr:::distinct(gene) %>%
						dplyr:::mutate(
							sample="s0",
							value = 0,
							ct = "bg",
							variable = "background"
						)
					else . 
				} %>%
				dplyr:::mutate_if(is.character, as.factor) %>%
				dplyr:::select(gene, ct, value) %>%
				dplyr:::group_by(gene, ct) %>%
				dplyr:::summarise(value = mean(value)) %>%
				dplyr:::ungroup() %>%
				tidyr:::spread(gene, value) %>%
				select(-ct)
		) %>%
		tibble::as_tibble() %>%
		dplyr:::mutate(sample = levels(beta_ancestor_run$sample)) %>%
		dplyr:::select(sample, everything())
		
	# balance the evidences
	e_map_matrix_stats = 
		e.obj$e %>% 
			dplyr:::group_by(ct_num) %>%
			dplyr:::summarise(dim = n()) %>%
			dplyr:::mutate(norm_factor = median(dim)/dim)
	# e_i_matrix still to solve
	
	my_local_design = if(is.null(my_design)) matrix(rep(1, ncol(mix)), ncol=1) else my_design

	# Create input object for the model
	model.in = list(
		G = fg %>% 
			dplyr:::distinct(gene) %>% 
			nrow(),
		
		S = mix %>% 
			dplyr:::distinct(sample) %>% 
			nrow(),
		
		P = fg %>% 
			dplyr:::distinct(ct) %>% 
			nrow(),
		
		R = my_local_design %>% 
			ncol(),
		
		y = mix %>% 
			dplyr:::select(gene, sample, value) %>% 
			tidyr:::spread(gene, value) %>%
			dplyr:::select(-sample), 
		
		X = my_local_design,
		
		x = fg %>% 
			dplyr:::select(gene, ct, value) %>%
			dplyr:::group_by(gene, ct) %>%
			dplyr:::summarise(value = mean(value)) %>%
			dplyr:::ungroup() %>%
			tidyr:::spread(ct, value) %>%
			dplyr:::select(-gene) %>%
			dplyr:::mutate_all(funs(ifelse(is.na(.), 0, .))),
		
		y_hat_background = y_hat_background %>%
			dplyr:::select(-sample),
		
		p_target = beta_ancestor_run %>%
			dplyr:::filter(ct==!!ct) %>%
			dplyr:::pull(value),
		
		theta = mix %>%
			dplyr:::distinct(sample, theta) %>%
			dplyr:::pull(theta),
		
		sigma_hyper_sd = sigma_hyper_sd,
		phi_hyper_sd = phi_hyper_sd,
		alpha_hyper_value = alpha_hyper_value,
		
		
		is_mix_microarray = as.numeric(is_mix_microarray),
		
		# For full Bayesian
		
		# Main node
		E =                               nrow(e.obj$e),
		E_MU =                            nrow(e.obj$e_mu),
		map_e_genes =                     e.obj$e$gene_num,
		map_e_ct =                        e.obj$e$ct_num,
		map_e_to_mu =                     e.obj$e$map_to_mu,
		e_ =                              e.obj$e$value,
		map_e_mu_gene =                   e.obj$e_mu$gene_num,
		map_e_mu_ct =                     e.obj$e_mu$ct_num,
		
		# e_i_matrix =                      e_i_matrix,
		# e_map_matrix_dim =                e_map_matrix_dim,
		# e_map_matrix_norm =               e_map_matrix_norm,
		
		
		# Background
		B =                               nrow(bg.obj$e),
		B_MU =                            nrow(bg.obj$e_mu),
		Q =                               length(unique(bg.obj$e$ct_num)),
		map_bg_genes =                    bg.obj$e$gene_num,
		map_bg_ct =                       bg.obj$e$ct_num,
		map_bg_to_mu =                    bg.obj$e$map_to_mu,
		b =                               bg.obj$e$value,
		map_bg_mu_gene =                  bg.obj$e_mu$gene_num,
		map_bg_mu_ct =                    bg.obj$e_mu$ct_num
	)

	browser()
	
	if(save_report) save(model.in, file=sprintf("%s/%s_model_in.RData", output_dir, ct))

	# Choose model
	model = if(fully_bayesian) stanmodels$ARMET_tc_recursive else stanmodels$ARMET_tcFix_recursive

	# Run model
	fit = 
		rstan::sampling(
			model,
			data=                             model.in,
			iter=                             1000 ,
			control =                         list(adapt_delta = 0.99, stepsize = 0.01, max_treedepth =15),
			cores=4
		)
	
	
	# Parse results

	proportions =                       t(parse_summary_vector_in_2D(apply( as.matrix(fit, pars = "beta"), 2, mean)))
	proportions_sd =                    t(parse_summary_vector_in_2D(apply( as.matrix(fit, pars = "beta"), 2, sd)))
	proportions_2.5 =                   t(parse_summary_vector_in_2D(apply( as.matrix(fit, pars = "beta"), 2, quantile, 0.025)))
	proportions_97.5 =                  t(parse_summary_vector_in_2D(apply( as.matrix(fit, pars = "beta"), 2, quantile, 0.975)))
	rownames(proportions) = rownames(proportions_2.5) = rownames(proportions_97.5) = rownames(proportions_sd) =  colnames(mix)
	colnames(proportions) =  colnames(proportions_2.5) = colnames(proportions_97.5) = colnames(proportions_sd) = order_cell_types

	# Save output
	if(save_report) save(fit, file=sprintf("%s/%s_fit.RData", output_dir, ct))
	p = rstan:::traceplot(fit, pars=c("alpha"), inc_warmup=F)
		
	#if(save_report) ggplot2:::ggsave(sprintf("%s_chains.png", ct), p)
	if(save_report)	write.csv(proportions, sprintf("%s/%s_composition.csv", output_dir, ct))

	
	list(proportions = proportions)
}
