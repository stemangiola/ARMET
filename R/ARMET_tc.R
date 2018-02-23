
#/
#|--------------------------------------------------------------------------
#| ARMET
#|--------------------------------------------------------------------------
#|
#| This function calls the stan model
#| Inpu: reference matrix, mix matrix.
#| Rownames are gene symbols, and colnames of ref are cell types
#\

# Initialize pipe
`%>%` <- magrittr::`%>%`

ARMET_tc = function(
	mix,
	my_design =                         NULL,
	cov_to_test =                       NULL,
	fully_bayesian =                    F,
	is_mix_microarray =                 F,
	observed_prop =                     NULL,
	ct_to_omit =                        c("t_CD4_naive", "adipocyte"),
	verbose =                           F,
	sigma_hyper_sd =                    0.1,
	phi_hyper_sd =                      2,
	alpha_hyper_value =                 2,
	save_report =                       F,
	custom_ref =                        NULL,
	multithread =                       T,
	do_debug =                          F,
	cell_type_root =                    "TME",
	choose_internal_ref =               NULL
){

	writeLines("ARMET: Started data processing")
	
	#Read ini file for some options
	get_ini()
	
	# Check input
	check_input(
		mix, 
		is_mix_microarray, 
		my_design, 
		cov_to_test, 
		sigma_hyper_sd, 
		custom_ref, 
		drop_node_from_tree(node_from_name(tree, cell_type_root), ct_to_omit) 
	)
	
	# Create directory
	output_dir = if(save_report)  create_temp_result_directory() else NULL

	# Format input
	mix = mix %>%	
		tidyr::gather("sample", "value", 2:ncol(mix)) %>%
		dplyr::mutate(gene = factor(gene), sample = factor(sample)) %>%
		{ if(max((.)$value) < 50) dplyr::mutate(value=exp(value)) else .}

	# Check if design matrix exists
	if(is.null(my_design)) 
		my_design = tibble::tibble(
			sample = levels(mix$sample),
			`(intercept)` = 1
		) %>%
		dplyr::mutate_if(is.character, as.factor)
		
	# Format tree
	my_tree =  format_tree( node_from_name(tree, cell_type_root), mix, ct_to_omit)

	# Ref formatting
	ref =
		{ 
			if(!is.null(custom_ref)) custom_ref
			else 
				if(!is_mix_microarray | (!is.null(choose_internal_ref) && choose_internal_ref == "ARNA")  ) ref_seq 
				else ref_array
		} %>% 
		tidyr::drop_na() %>% 
		dplyr::filter(
			ct %in% 
			get_leave_label(my_tree, last_level = 0, label = "name")
		) %>% 
		droplevels()
	
	
	# Calculate stats for ref
	if(save_report) write.csv(get_stats_on_ref(ref, my_tree), sprintf("%s/stats_on_ref.csv", output_dir))

	# Make data sets comparable
	common_genes =                      intersect(as.character(mix$gene), as.character(ref$gene))
	mix =                               mix %>% dplyr::filter(gene%in%common_genes) %>% droplevels()
	ref =                               ref %>% dplyr::filter(gene%in%common_genes) %>% droplevels()

	# Normalize data
	norm.obj = 													wrapper_normalize_mix_ref(mix, ref, is_mix_microarray)
	ref = 															norm.obj$ref 
	mix = 															norm.obj$mix
	
	# Round if RNA seq
	if(!is_mix_microarray) ref = 				ref %>% dplyr::mutate(value=round(value))
	if(!is_mix_microarray) mix = 				mix %>% dplyr::mutate(value=round(value))
	
	# Plot densities
	#plot_densities(df)
	
	# Save density plot
	if(save_report) ggplot2::ggsave(sprintf("%s/densities.png", output_dir), plot=norm.obj$plot)
	
	## Execute core ##############################################################################
	##############################################################################################
	
	my_tree = 
		run_coreAlg_though_tree(
			my_tree, 
			list(
				mix =                           mix, 
				ref =                           ref, 
				my_design=                      my_design, 
				cov_to_test =                   cov_to_test, 
				fully_bayesian =                fully_bayesian,
				observed_prop =                 observed_prop, 
				ct_to_omit =                    ct_to_omit, 
				my_tree =                       my_tree, 
				is_mix_microarray =             is_mix_microarray,
				save_report =                   save_report,
				output_dir =                    output_dir,
				sigma_hyper_sd =                sigma_hyper_sd,
				phi_hyper_sd =                  phi_hyper_sd,
				alpha_hyper_value =             alpha_hyper_value,
				multithread =                   multithread,
				do_debug =                      do_debug
			)
	)

	##############################################################################################
	##############################################################################################

	# Create tree with hypothesis testing
	osNode.stat = 
		switch(
			(!is.null(cov_to_test)) + 1,
			NULL,
			get_tree_hypoth_test(osNode, my_tree)
		)
		
	if(save_report) save(osNode.stat, file=sprintf("%s/tree_pvalues.RData", output_dir))

	# Return
	list(
		
		proportions =	get_last_existing_leaves_with_annotation( my_tree ) %>%
			dplyr::select(-relative_proportion) %>%
			tidyr::spread(ct, absolute_proportion),
		
		signatures = 
			list(
				orig =  
					ref %>% 
					dplyr::select(gene, ct, value) %>%
					dplyr::group_by(gene, ct) %>%
					dplyr::summarise(value = mean(value)) %>%
					dplyr::ungroup() %>%
					tidyr::spread(ct, value),
				predicted = NULL
			),
		
		mix = mix %>% 
			dplyr::filter(gene %in% get_genes( my_tree )) %>%
			droplevels(),
		
		stats = osNode.stat
	)
	
}
