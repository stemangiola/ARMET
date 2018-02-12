
#/
#|--------------------------------------------------------------------------
#| ARMET
#|--------------------------------------------------------------------------
#|
#| This function calls the stan model
#| Inpu: reference matrix, mix matrix.
#| Rownames are gene symbols, and colnames of ref are cell types
#\

ARMET_tc = function(
	mix,
	my_design =                         NULL,
	cov_to_test =                       NULL,
	fully_bayesian =                    F,
	is_mix_microarray =                 F,
	observed_prop =                     NULL,
	ct_to_omit =                        c("t_CD4_naive", "adipocyte"),
	verbose =                           F,
	sigma_hyper_sd =                    0.01,
	phi_hyper_sd =                      5,
	alpha_hyper_value =                 10,
	save_report =                       F,
	custom_ref =                        NULL
){
	
	#Read ini file for some options
	get_ini()
	
	# Check input
	check_input(mix, is_mix_microarray, my_design, cov_to_test, sigma_hyper_sd, custom_ref, tree)
	
	# Create directory
	output_dir = if(save_report)        create_temp_result_directory() else NULL
	
	# Format input
	mix =                               as.matrix(mix)
	if(max(mix) < 50) mix =             exp(mix)
	
	##########################################################

	# Load reference
	# ref_RNAseq = as.matrix(read.csv("~/PhD/deconvolution/ARMET_dev/ARMET_TME_signature_df_RNAseq.csv", header=T, row.names=1))
	# colnames(ref_RNAseq) = as.vector(sapply(colnames(ref_RNAseq), function(cn) strsplit(cn, ".", fixed=T)[[1]][1]))
	# save(ref_RNAseq, file="data/ref_RNAseq.rda")
	# ref_array = as.matrix(read.csv("~/PhD/deconvolution/ARMET_dev/ARMET_TME_signature_df_array.csv", header=T, row.names=1))
	# colnames(ref_array) = as.vector(sapply(colnames(ref_array), function(cn) strsplit(cn, ".", fixed=T)[[1]][1]))
	# save(ref_array, file="data/ref_array.rda")
	if(!is.null(custom_ref))            ref = custom_ref 
	else                                ref = if(!is_mix_microarray) ref_RNAseq else ref_array

	
	# Create trees
	# library(jsonlite)
	# tree = read_json("/wehisan/home/allstaff/m/mangiola.s/PhD/deconvolution/ARMET_dev/ARMET_TME_tree_RNAseq.json")
	# save(tree, file="data/tree_json.rda")
	data(tree_json)
	my_tree =                           drop_node_from_tree(tree, ct_to_omit)
	ref =                               ref[, colnames(ref)%in%get_leave_label(my_tree, last_level = 0, label = "name")]
	
	# Make data sets comparable
	common_genes =                      intersect(rownames(mix), rownames(ref))
	mix =                               mix[common_genes,, drop=FALSE]
	ref =                               ref[common_genes,, drop=FALSE]
	
	# Normalize data
	norm.obj = 													wrapper_normalize_mix_ref(mix, ref, is_mix_microarray)
	ref = 															norm.obj$ref
	mix = 															norm.obj$mix
	
	# Save density plot
	if(save_report) ggplot2:::ggsave(sprintf("%s/densities.png", output_dir), plot=norm.obj$plot)
	
	## Execute core ##############################################################################
	##############################################################################################
	
	my_tree = run_coreAlg_though_tree(
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
			alpha_hyper_value =             alpha_hyper_value
		)
	)
		
	##############################################################################################
	##############################################################################################

	# Create one tree per sample
	my_trees =                            divide_trees_proportion_across_many_trees(my_tree)
	
	# Calculate absolute proportions
	my_trees = 														add_absolute_proportions_to_trees(my_trees)
	
	# Produce table proportions
	proportions = 											  get_proportions_table(my_trees)
	
	# Check if markers were reliable
	markers =                             intersect(intersect(get_genes(tree), rownames(ref)), rownames(mix))
	signatures = 												  list()
	signatures$orig = 									  get_mean_signature(ref[markers,], verbose = F)
	signatures$orig =                     signatures$orig[,colnames(signatures$orig)%in%colnames(proportions)]
	signatures$predicted =							  NULL
	
	# Create tree with hypothesis testing
	if(!is.null(cov_to_test) & 0) {
		data(osNode)
		tree.test = get_tree_hypoth_test(osNode, path=output_dir, cov_to_test)
		print(tree.test, "significance_causal", "pvalue_causal", "direction_causal", "significance_effectual", "pvalue_effectual", "direction_effectual")
		save(tree.test, file=sprintf("%s/tree_pvalues.RData", output_dir))
	}
	
	return(list(
		proportions =                       proportions,
		signatures =                        signatures,
		mix =                               mix[markers,,drop=F]
	))
	
}
