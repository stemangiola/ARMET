#' Core algorithm
#'
#' @description Calls the Stan model
#'
#' @importFrom rstan sampling
#' @importFrom rstan summary
#'
#' @importFrom tibble tibble
#' @importFrom tibble as_tibble
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom dplyr rename
#' @importFrom dplyr distinct
#' @importFrom dplyr filter
#' @importFrom dplyr mutate_if
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom dplyr summarise
#' @importFrom dplyr pull
#' @importFrom dplyr bind_rows
#' @importFrom dplyr mutate_all
#' @importFrom dplyr funs
#' @importFrom dplyr everything
#' @importFrom dplyr arrange
#' @importFrom dplyr bind_cols
#' @importFrom dplyr one_of
#'
#' @importFrom tidyr spread
#' @importFrom tidyr gather
#'
#' @importFrom tidybayes gather_samples
#' @importFrom tidybayes spread_samples
#' @importFrom tidybayes median_qi
#' @importFrom tidybayes mean_qi
#'
#' @importFrom abind abind
#'
#'
#'
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
	save_fit =              obj.in$save_fit
	seed =                  obj.in$seed
	phi =                   obj.in$phi

	# Get ref of the current level
	ref =
		ref %>%
		# Get cell type just for this recursive peers and descendants
		mutate(ct = as.character(ct)) %>%
		left_join(
			get_map_foreground_background(tree, ct) %>%
				select(-ct) %>%
				rename(ct=ancestor) %>%
				distinct(),
			by="ct"
		) %>%
		# Filter non used cell types
		filter(!is.na(variable)) %>%
		# Filter non used genes
		filter(
			gene %in%
			get_node_label_level_specfic(
				get_node_from_name(my_tree,  !!ct),
				label = "markers",
				start_level = 1,
				stop_level = 1
			)
		) %>%
		mutate_if(is.character, as.factor) %>%
		droplevels()

	# Print stats on ref
	get_stats_on_ref(ref,my_tree) %>%
		filter(ct %in% unique(ref$ct))

	# filter mix and add theta value
	mix = mix %>%
		group_by(sample) %>%
		mutate(
			theta=
				length(which(value>0)) /
				n()
		) %>%
		ungroup() %>%
		filter(gene %in% unique(ref$gene))

	fg = ref %>%
		filter(variable=="main") %>%
		droplevels()

	# Setup background
	bg = 	ref %>%
		filter(variable=="background") %>%
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
			group_by(sample) %>%
			summarise(tot=sum(absolute_proportion)) %>%
			pull(tot) != 1
		)
	) {
		writeLines("ARMET: The absolute proportions are supposed to sum to 1 for each sample")
		print(
			ancestor_run_prop_table %>%
				group_by(sample) %>%
				summarise(tot=sum(absolute_proportion)) %>%
				filter(tot != 1)
		)
		#stop()
	}

	# Get the probability table of my cell type
	fg_prop = ancestor_run_prop_table %>%
		filter(ct==!!ct) %>%
		droplevels()

	# Get the probability table of the background
	bg_prop = ancestor_run_prop_table %>%
		filter(ct != !!ct) %>%
		droplevels()

	# Prepare input for full bayesian
	# e.obj = prepare_input(fg, my_tree)
	# bg.obj = prepare_input(bg, my_tree)

	# Calculate the value of the genes for background
	# !!rlang::sym(ct) -> for variable name without quotes

	y_hat_background =
		as.matrix(
			bg_prop %>%
				switch(
					(!nrow(.)==0) + 1,
					(.) %>% bind_rows(
						fg_prop %>%
							mutate(
								ct = factor("bg"),
								relative_proportion = 0,
								absolute_proportion = 0
							)
					),
					.
				) %>%
				select(sample, ct, absolute_proportion) %>%
				spread(ct, absolute_proportion) %>%
				select(-sample)
		) %*%
		as.matrix(
			bg %>%
			{
				if(nrow(.)==0)
					# Take fg as template
					fg %>%
					distinct(gene) %>%
					mutate(
						sample="s0",
						value = 0,
						ct = "bg",
						variable = "background"
					)
				else .
			} %>%
				mutate_if(is.character, as.factor) %>%
				select(gene, ct, value) %>%
				spread(gene, value) %>%
				select(-ct) %>%
				mutate_all(funs(ifelse(is.na(.), 0, .)))
		) %>%
		as_tibble() %>%
		mutate(sample = levels(mix$sample)) %>%
		select(sample, everything())

	# Create input object for the model
	model.in = list(

		G = fg %>%
			distinct(gene) %>%
			nrow(),

		S = mix %>%
			distinct(sample) %>%
			nrow(),

		P = fg %>%
			distinct(ct) %>%
			nrow(),

		R = my_design %>%
			select(-sample) %>%
			ncol(),

		y = mix %>%
			select(gene, sample, value) %>%
			spread(gene, value) %>%
			arrange(sample) %>%
			select(-sample),

		X = my_design %>%
			arrange(sample) %>%
			select(-sample) %>%
			# transform factors into numeric safely
			mutate_if(is.factor, as.character) %>%
			mutate_if(is.character, as.numeric) %>%
			mutate_if(any_column_double, scale),

		x_genes = fg %>%
			pull(gene) %>%
			levels() %>%
			as_tibble() %>%
			rename(gene = value),

		x = fg %>%
			select(gene, ct, value) %>%
			spread(ct, value) %>%
			mutate_if(is.numeric, funs(ifelse(is.na(.), 0, .))) %>%
			select(-gene),

		y_hat_background = y_hat_background %>%
			select(-sample),

		p_ancestors = ancestor_run_prop_table %>%
			select(sample, ct, absolute_proportion) %>%
			spread(ct, absolute_proportion) %>%
			select(-sample) %>%
			select(-one_of(ct), one_of(ct)),

		P_a = length(levels(ancestor_run_prop_table$ct)),

		theta = mix %>%
			distinct(sample, theta) %>%
			pull(theta) %>%
			as.array(),

		is_mix_microarray =  as.numeric(is_mix_microarray),
		omit_regression =    as.numeric(omit_regression),

		# Horseshoe
		nu_local = 1,
		nu_global = 1,
		par_ratio = 0.8,
		slab_df = 4,
		slab_scale = 0.5

	)

	# Save the raw model resul for debugging
	if(save_report) save(model.in, file=sprintf("%s/%s_model_in.RData", output_dir, ct))

	# Choose model
	model = switch(is_mix_microarray + 1, stanmodels$ARMET_tcFix_recursive, stanmodels$ARMET_tcFix_recursive_array)

	# Run model
	fit =
		rstan::sampling(
			model,
			data=                             model.in,
			iter=                             ifelse(ct %in% c("TME", "immune_cell"), 600, 1400) ,
			warmup =                          ifelse(ct %in% c("TME", "immune_cell"), 400, 700),
			#control =                         list(adapt_delta = 0.9, stepsize = 0.01, max_treedepth =15),
			#control =                         list(max_treedepth =15),
			cores = 4,
			seed = ifelse(is.null(seed), sample.int(.Machine$integer.max, 1), seed),
			save_warmup=FALSE
		)

	# Add names ct in this node
	node =
		add_info_to_tree(
			node,
			ct,
			"ct_in_analysis",
			c(colnames(model.in$p_ancestors)[-ncol(model.in$p_ancestors)], colnames(model.in$x))
		)

	# Parse results
	proportions =
		fit %>%
		gather_samples(beta[sample_idx, ct_idx]) %>%
		mean_qi() %>%
		ungroup() %>%
		annotate_posterior_sample_ct(colnames(model.in$x),my_design) %>%
		rename(relative_proportion = estimate) %>%
		mutate_if(is.character, as.factor)

	# Add info to the node
	node = add_data_to_tree_from_table(node, proportions, "relative_proportion", append = T)

	# Set up background trees
	node = add_absolute_proportions_to_tree(node)

	# Add hypothesis testing
	node = switch(
		(!save_fit) + 1,
		add_info_to_tree(	node,	ct,	"fit",fit,	append = F),
		node
	)

	node = switch(
		is.null(cov_to_test) + 1,
		add_info_to_tree(
			node,
			ct,
			"estimate_prop_with_uncertanties",
			fit %>%
				gather_samples(beta_global[sample_idx, ct_idx]) %>%
				median_qi() %>%
				annotate_posterior_sample_ct(node$ct_in_analysis,my_design),
			append = F
		),
		node
	)

	node = switch(
		is.null(cov_to_test) + 1,
		add_info_to_tree(
			node,
			ct,
			"estimate_generated_prop_with_uncertanties",
			fit %>%
				gather_samples(beta_gen[sample_idx, ct_idx]) %>%
				median_qi() %>%
				annotate_posterior_sample_ct(node$ct_in_analysis,my_design),
			append = F
		),
		node
	)

	# Beta regression
	fit_betaReg =
		rstan::sampling(
			stanmodels$ARMET_betaReg, # rstan::stan_model("src/stan_files/ARMET_betaReg.stan"), #
			data= list(
				S = model.in$S,
				P = model.in$P,
				R = model.in$R,
				X = model.in$X,
				N = 100,

				# Take posterior to a 3D matrix
				beta = (
					fit %>%
						gather_samples(beta_global[sample_idx, ct_idx]) %>%
						filter(ct_idx > ncol(model.in$p_ancestors)-1) %>%
						filter(.iteration <= 25) %>%
						select(estimate, sample_idx, ct_idx, .iteration) %>%
						mutate(.iteration = 1:n()) %>%
						ungroup() %>%
						select(-term) %>%

						# Adjust for lower and upper limits
						#mutate(estimate = (estimate*(100-1) + (1/model.in$P) ) / (100)) %>%

						# Convert to 3D
						{
							foreach(i = (.) %>% pull(ct_idx) %>% unique()) %do%
							{
								(.) %>% filter(ct_idx==i) %>% select(-ct_idx) %>% spread( sample_idx, estimate) %>% select(-.iteration)
							}
						} %>%
						abind(along=3) %>%
						aperm(c(2,3,1))
				)
			),
			cores = 4,
			seed = ifelse(is.null(seed), sample.int(.Machine$integer.max, 1), seed),
			save_warmup=FALSE,
			iter=     ifelse(ct %in% c("TME", "immune_cell"), 600, 1400) ,
			warmup =  ifelse(ct %in% c("TME", "immune_cell"), 400, 700)
		)


	# Hypothesis test
	node = switch(
		is.null(cov_to_test) + 1,
		add_info_to_tree(
			node,
			ct,
			"stats",
			left_join(

				# Extrinsic regression
				fit %>%
					spread_samples(extrinsic[covariate_idx, ct_idx]) %>%
					median_qi() %>%
					ungroup() %>%
					left_join(
						tibble(
							covariate = colnames(my_design %>% select(-sample)),
							covariate_idx = 1:ncol(my_design %>% select(-sample))
						),
						by = "covariate_idx"
					) %>%
					left_join(
						tibble(
							ct = node$ct_in_analysis,
							ct_idx = 1:length(node$ct_in_analysis)
						),
						by = "ct_idx"
					) %>%
					filter(ct %in% colnames(model.in$x)) %>%
					mutate(e.sig = ifelse(conf.low * conf.high > 0, "*", "")) %>%
					setNames(gsub("conf", "e", colnames(.))) %>%
					select(-covariate_idx, -ct_idx, -.prob),

				# Extrinsic regression
				fit_betaReg %>%
					spread_samples(intrinsic[covariate_idx, ct_idx]) %>%
					median_qi() %>%
					ungroup() %>%
					left_join(
						tibble(
							covariate = colnames(my_design %>% select(-sample)),
							covariate_idx = 1:ncol(my_design %>% select(-sample))
						),
						by = "covariate_idx"
					) %>%
					left_join(
						tibble(
							ct = colnames(model.in$x),
							ct_idx = 1:ncol(model.in$x)
						),
						by = "ct_idx"
					) %>%
					mutate(i.sig = ifelse(conf.low * conf.high > 0, "*", "")) %>%
					setNames(gsub("conf", "i", colnames(.))) %>%
					select(-covariate_idx, -ct_idx, -.prob),
				by = c("covariate",  "ct")
			) %>% select(covariate, ct, everything())
		),
		node
	)

	# Save output
	if(save_report) save(fit, file=sprintf("%s/%s_fit.RData", output_dir, ct))
	if(save_report)	write.csv(proportions, sprintf("%s/%s_composition.csv", output_dir, ct))


	list(proportions = proportions, node = node)
}
