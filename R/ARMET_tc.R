#' ARMET-tc main
#'
#' @description This function calls the stan model.
#'
#'
#' @importFrom tibble tibble
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom dplyr mutate_if
#'
#' @importFrom tidyr spread
#' @importFrom tidyr gather
#' @importFrom tidyr drop_na
#'
#' @importFrom tidybayes gather_samples
#' @importFrom tidybayes median_qi
#'
#' @importFrom foreach foreach
#' @importFrom foreach %do%
#'
#' @param mix A matrix
#' @param my_design A matrix
#' @param cov_to_test A character string
#' @param fully_bayesian A boolean
#' @param is_mix_microarray A boolean
#' @param ct_to_omit A character string
#' @param verbose A boolean
#' @param save_report A boolean
#' @param custom_ref A matrix
#' @param multithread A boolean
#' @param do_debug A boolean
#' @param cell_type_root A character string
#' @param choose_internal_ref A design matrix
#' @param omit_regression A boolean
#' @param save_fit A boolean
#' @param seed An integer
#'
#' @return An ARMET object
#'
#' @export
#'
ARMET_tc = function(
	mix,
	my_design =                         NULL,
	cov_to_test =                       NULL,
	is_mix_microarray =                 F,
	ct_to_omit =                        c("t_CD4_naive", "adipocyte"),
	verbose =                           F,
	save_report =                       F,
	custom_ref =                        NULL,
	multithread =                       T,
	do_debug =                          F,
	cell_type_root =                    "TME",
	choose_internal_ref =               NULL,
	omit_regression =                   F,
	save_fit =                          F,
	seed =                              NULL
){

	input = c(as.list(environment()))

	writeLines("ARMET: Started data processing")

	#Read ini file for some options
	get_ini()

	# Check input
	check_input(
		mix,
		is_mix_microarray,
		my_design,
		cov_to_test,
		custom_ref,
		drop_node_from_tree(get_node_from_name(tree, cell_type_root), ct_to_omit)
	)

	# Create directory
	output_dir = if(save_report)  create_temp_result_directory() else NULL

	# Format input
	mix = mix %>%
		gather("sample", "value", 2:ncol(mix)) %>%
		mutate_if(is.factor, as.character) %>%
		mutate_if(is.character, as.factor) %>%
		mutate_if(is.factor, factor) %>%
		{ if(max((.)$value) < 50) mutate(value=exp(value)) else .}

	# Check if design matrix exists
	my_design =
		switch(
			is.null(my_design) + 1,
			my_design,
			tibble(	sample = levels(mix$sample),	`(intercept)` = 1	)
		) %>%
		mutate_if(is.factor, as.character) %>%
		mutate_if(is.character, as.factor) %>%
		mutate_if(is.factor, factor)

	# Format tree
	my_tree =  format_tree( get_node_from_name(tree, cell_type_root), mix, ct_to_omit)

	# Ref formatting
	ref =
		# See if custom ref is provided
		switch(
			is.null(custom_ref) + 1,
			ref_to_summary_ref(my_tree, custom_ref),
			# See if choose internal ref is provided
			switch(
				is.null(choose_internal_ref) + 1,
				switch(
					(choose_internal_ref == "ARNA") + 1,
					ref_array_recursive,
					ref_RNAseq_recursive
				),
				# if npthing set choose default
				switch(
					(!is_mix_microarray) + 1,
					ref_array_recursive,
					ref_RNAseq_recursive
				)
			)
		) %>%
		drop_na() %>%
		filter(
			ct %in%
			get_leave_label(my_tree, last_level = 0, label = "name")
		) %>%
		droplevels()


	# Calculate stats for ref
	if(save_report) write.csv(get_stats_on_ref(ref, my_tree), sprintf("%s/stats_on_ref.csv", output_dir))

	# Make data sets comparable
	common_genes =                      intersect(as.character(mix$gene), as.character(ref$gene))
	mix =                               mix %>% filter(gene%in%common_genes) %>% droplevels()
	ref =                               ref %>% filter(gene%in%common_genes) %>% droplevels()

	# Normalize data
	mix = 													wrapper_normalize_mix_ref(mix, ref, is_mix_microarray)

	# Round if RNA seq
	if(!is_mix_microarray) ref = 				ref %>% mutate(value=round(value))
	if(!is_mix_microarray) mix = 				mix %>% mutate(value=round(value))

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
				ct_to_omit =                    ct_to_omit,
				my_tree =                       my_tree,
				is_mix_microarray =             is_mix_microarray,
				save_report =                   save_report,
				output_dir =                    output_dir,
				multithread =                   multithread,
				do_debug =                      do_debug,
				omit_regression =               omit_regression,
				save_fit =                      save_fit,
				seed =                          seed,
				verbose =                       verbose
			)
	)

	##############################################################################################
	##############################################################################################

	writeLines("ARMET: building output")

	# # Create tree with hypothesis testing
	# library(data.tree)
	# fileName=("~/PhD/deconvolution/ARMET_dev/ARMET_TME_tree.yaml")
	# ya=readChar(fileName, file.info(fileName)$size)
	# osList <- yaml::yaml.load(ya)
	# treeYaml <- as.Node(osList)
	# save(treeYaml, file="data/treeYaml.rda")

	treeYaml.stat =
		switch(
			(!is.null(cov_to_test)) + 1,
			NULL,
			get_tree_hypoth_test(treeYaml, my_tree, cov_to_test)
		)

	if(save_report) save(treeYaml.stat, file=sprintf("%s/tree_pvalues.RData", output_dir))

	# Return
	list(

		# Matrix of proportions
		proportions =	get_last_existing_leaves_with_annotation( my_tree ) %>%
			select(-relative_proportion) %>%
			spread(ct, absolute_proportion),

		# What mixture was used by the model after normalization
		mix = mix %>%
			filter(gene %in% get_genes( my_tree )) %>%
			droplevels(),

		# Return the statistics
		stats = treeYaml.stat,

		# Return the annotated tree
		tree = my_tree,

		# Return the input itself
		input = input
	)

}
