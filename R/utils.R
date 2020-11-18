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
		axis.text.x = element_text(
			angle = 90,
			hjust = 1,
			vjust = 0.5
		),
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
		))
	)

# Greater than
gt = function(a, b){	a > b }

# Smaller than
st = function(a, b){	a < b }

# Negation
not = function(is){	!is }

#' This is a generalisation of ifelse that acceots an object and return an objects
#' @keywords internal
#'
#' @import dplyr
#' @import tidyr
#'
#' @param input.df A tibble
#' @param condition A boolean
#' @return A tibble
ifelse_pipe = function(.x, .p, .f1, .f2 = NULL) {
	switch(.p %>% `!` %>% sum(1),
				 as_mapper(.f1)(.x),
				 if (.f2 %>% is.null %>% `!`)
				 	as_mapper(.f2)(.x)
				 else
				 	.x)

}

#' This is a generalisation of ifelse that acceots an object and return an objects
#' @keywords internal
#'
#' @import dplyr
#' @import tidyr
#'
#' @param .x A tibble
#' @param .p1 A boolean
#' @param .p2 ELSE IF condition
#' @param .f1 A function
#' @param .f2 A function
#' @param .f3 A function
#'
#' @return A tibble
ifelse2_pipe = function(.x, .p1, .p2, .f1, .f2, .f3 = NULL) {
	# Nested switch
	switch(# First condition
		.p1 %>% `!` %>% sum(1),
		
		# First outcome
		as_mapper(.f1)(.x),
		switch(
			# Second condition
			.p2 %>% `!` %>% sum(1),
			
			# Second outcome
			as_mapper(.f2)(.x),
			
			# Third outcome - if there is not .f3 just return the original data frame
			if (.f3 %>% is.null %>% `!`)
				as_mapper(.f3)(.x)
			else
				.x
		))
}

#' This is a generalisation of ifelse that acceots an object and return an objects
#' @keywords internal
#'
#' @import dplyr
#' @import tidyr
#'
#' @param .x A tibble
#' @param .p1 A boolean
#' @param .p2 ELSE IF condition
#' @param .f1 A function
#' @param .f2 A function
#' @param .f3 A function
#'
#' @return A tibble
ifelse3_pipe = function(.x, .p1, .p2, .p3, .f1, .f2, .f3, .f4 = NULL) {
	# Nested switch
	switch(# First condition
		.p1 %>% `!` %>% sum(1),
		
		# First outcome
		as_mapper(.f1)(.x),
		switch(
			# Second condition
			.p2 %>% `!` %>% sum(1),
			
			# Second outcome
			as_mapper(.f2)(.x),
			
			# Third outcome - if there is not .f3 just return the original data frame
			switch(
				# Second condition
				.p3 %>% `!` %>% sum(1),
				
				# Second outcome
				as_mapper(.f3)(.x),
				
				# Third outcome - if there is not .f3 just return the original data frame
				if (.f4 %>% is.null %>% `!`)
					as_mapper(.f4)(.x)
				else
					.x
			)
		))
}

#' This is a generalisation of ifelse that acceots an object and return an objects
#' @keywords internal
#'
#' @import dplyr
#' @import tidyr
#'
#' @param .x A tibble
#' @param .p1 A boolean
#' @param .p2 ELSE IF condition
#' @param .f1 A function
#' @param .f2 A function
#' @param .f3 A function
#'
#' @return A tibble
ifelse4_pipe = function(.x, .p1, .p2, .p3, .p4, .f1, .f2, .f3, .f4, .f5 = NULL) {
	# Nested switch
	switch(# First condition
		.p1 %>% `!` %>% sum(1),
		
		# First outcome
		as_mapper(.f1)(.x),
		switch(
			# Second condition
			.p2 %>% `!` %>% sum(1),
			
			# Second outcome
			as_mapper(.f2)(.x),
			
			# Third outcome - if there is not .f3 just return the original data frame
			switch(
				# Second condition
				.p3 %>% `!` %>% sum(1),
				
				# Second outcome
				as_mapper(.f3)(.x),
				
				# Third outcome - if there is not .f3 just return the original data frame
				switch(
					# Second condition
					.p4 %>% `!` %>% sum(1),
					
					# Second outcome
					as_mapper(.f4)(.x),
					
					# Third outcome - if there is not .f3 just return the original data frame
					if (.f5 %>% is.null %>% `!`)
						as_mapper(.f5)(.x)
					else
						.x
				)
			)
		))
}

# @description Add partition column dto data frame
add_partition = function(df.input, partition_by, n_partitions) {
	df.input %>%
		left_join(
			(.) %>%
				select(!!partition_by) %>%
				distinct %>%
				mutate(
					partition = 1:n() %>%
						divide_by(length((.))) %>%
						#	multiply_by(min(n_partitions, df.input %>% distinct(symbol) %>% nrow)) %>%
						multiply_by(n_partitions) %>%
						ceiling
				)
		)
}

# @description  Plot differences between inferred and observed transcription abundances
plot_differences_in_lambda = function() {
	# Plot differences in lambda_log
	(
		fit %>%
			tidybayes::gather_draws(lambda_log[G]) %>%
			tidybayes::median_qi() %>%
			left_join(
				counts_baseline %>%
					distinct(symbol, G, `Cell type category`)
			) %>%
			left_join(
				reference_filtered %>%
					distinct(symbol, lambda_log, `Cell type category`) %>%
					rename(`lambda_log` = lambda_log)
			) %>%
			ggplot(aes(
				x = lambda_log, y = .value, label = G
			)) + geom_point() + geom_abline(
				intercept = 0,
				slope = 1,
				color = "red"
			) + my_theme
	)  
	#%>% plotly::ggplotly()

	#
	(
		fit %>%
			extract(
				pars = c(
					"lambda_mu",
					"lambda_sigma",
					# "exposure_rate",
					"lambda_log",
					"sigma_inv_log",
					"prop"
				)
			) %>%
			as.data.frame %>% as_tibble() %>%
			mutate(chain = rep(1:3, 100) %>% sort %>% as.factor) %>%
			select(chain, everything()) %>% gather(par, draw, -chain) %>%
			group_by(chain, par) %>%
			summarise(d = draw %>% median) %>%
			ggplot(aes(
				y = d, x = par, color = chain
			)) + geom_point()
	) 
	#%>% plotly::ggplotly()

}

# @description Print which reference genes are in the mix
get_overlap_descriptive_stats = function(mix_tbl, ref_tbl) {
	writeLines(
		sprintf(
			"%s house keeping genes are missing from the input mixture",
			ref_tbl %>% filter(`house keeping`) %>% distinct(symbol) %>% anti_join(mix_tbl %>% distinct(symbol)) %>% nrow
		)
	)

	writeLines(
		sprintf(
			"%s marker genes are missing from the input mixture",
			ref_tbl %>% filter(!`house keeping`) %>% distinct(symbol) %>% anti_join(mix_tbl %>% distinct(symbol)) %>% nrow
		)
	)

	ref_tbl %>%
		filter(!`house keeping`) %>% distinct(symbol, ct1, ct2) %>% anti_join(mix_tbl %>% distinct(symbol)) %>%
		count(ct1, ct2) %>%
		rename(`missing markers` = n) %>%
		print(n = 999)

}

#@description Get data format for MPI deconvolution part
plot_counts_inferred_sum = function(fit_obj, samples = NULL, level) {
	fit_obj$internals$fit[[level]] %>%
		rstan::summary(par = c("nb_sum")) %$% summary %>%
		as_tibble(rownames = "par") %>% select(par, `2.5%`, `50%`, `97.5%`) %>%
		separate(par, c(".variable", "Q", "GM"), sep = "\\[|,|\\]") %>%
		mutate(Q = Q %>% as.integer, GM = GM %>% as.integer) %>%
		left_join(fit_obj %$% data_source %>% distinct(Q, GM, symbol, count, sample)) %>%
 
		# Select samples
		{
			if (samples %>% is.null %>% `!`)
				(.) %>% filter(sample %in% !samples)
			else
				(.)
		} %>%

		# Check if inside
		rowwise %>%
		mutate(inside = between(count, `2.5%`, `97.5%`)) %>%
		ungroup %>%
		ggplot(aes(
			x = count + 1,
			y = `50%` + 1,
			color = inside,
			label = symbol
		)) +
		geom_point(alpha = 0.5)  +
		geom_abline(slope = 1, intercept = 0) +
		geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.5) +
		facet_wrap(~ sample) +
		scale_y_log10() +
		scale_x_log10()

}

#@description Get which chain cluster is more opulated in case I have divergence
choose_chains_majority_roule = function(fit_parsed) {
	fit_parsed %>%
		inner_join(
			# Calculate modes
			fit_parsed %>%
				select(.chain, .value) %>%
				{
					main_cluster =
						(.) %>%
						pull(.value) %>%
						kmeans(centers = 2) %>% {
							bind_cols(
								(.) %$% centers %>% as_tibble(.name_repair = "minimal") %>% setNames("center") ,
								(.) %$% size %>% as_tibble(.name_repair = "minimal")  %>% setNames("size")
							)
						} %>%
						arrange(size %>% desc) %>%
						slice(1)

					(.) %>%
						group_by(.chain) %>%
						summarise(
							.lower_chain = quantile(.value, probs = c(0.025)),
							.upper_chain = quantile(.value, probs = c(0.975))
						) %>%
						ungroup %>%
						mutate(center = main_cluster %>% pull(center))
				} %>%

				# Filter cains
				rowwise() %>%
				filter(between(center, .lower_chain, .upper_chain)) %>%
				ungroup %>%
				distinct(.chain),
			by = ".chain"

		)
}

#@description Filter the reference
filter_reference = function(reference, mix, n_markers) {

	# Check if all cell types in ref are in n_markers
	if( reference %>% filter(!ct1 %in% n_markers$ct1) %>% nrow %>% `>` (0) )
		stop(
			sprintf(
				"The cell types %s are present in reference but not in n_markers dataset",
				reference %>% filter(!ct1 %in% n_markers$ct1) %>% distinct(ct1) %>% pull(1) %>% paste(collapse=", ")
				)
		)

	reference %>%
		filter(symbol %in% (mix %>% distinct(symbol) %>% pull(symbol) )) %>%
		{
			bind_rows(
				# Get markers based on common genes with mix
				(.) %>%
					filter(!`house keeping`) %>%
					group_by(ct1, ct2) %>%
					do({
						n_markers = (.) %>% slice(1) %>% pull(`n markers`)
						(.) %>%
							inner_join(
								(.) %>%
									distinct(symbol, rank) %>%
									arrange(rank) %>%
									slice(1:n_markers),
								by = c("symbol", "rank")
							)

					}) %>%
					ungroup() %>%

					# Filter markers in upper levels that have been selected for lower levels
					inner_join(
						(.) %>%
							distinct(symbol, level) %>%
							group_by(symbol) %>%
							arrange(level %>% desc) %>%
							slice(1) %>%
							ungroup,
						by = c("level", "symbol")
					),

					# %>%
					#
					# # Print number of markewrs per comparison
					# {
					# 	(.) %>% distinct(symbol, ct1, ct2) %>% count(ct1, ct2) %>% print(n = 99)
					# 	(.)
					# },

				# Get house keeping genes
				(.) %>% filter(`house keeping`)
			)
		}
}

# get_idx_level
get_idx_level = function(tree, my_level) {
	left_join(
		tree %>% data.tree::ToDataFrameTree("name") %>% as_tibble,
		data.tree::Clone(tree) %>%
			{
				data.tree::Prune(., function(x)
					x$level <= my_level + 1)
				.
			} %>%
			data.tree::ToDataFrameTree("level", "C", "isLeaf", "name") %>%
			as_tibble %>%
			filter(isLeaf) %>%
			left_join(
				tree %>%
					data.tree::ToDataFrameTree("name", "level", "isLeaf") %>%
					as_tibble %>%
					filter(level <= my_level & isLeaf) %>%
					mutate(isAncestorLeaf = !isLeaf) %>%
					select(name, isAncestorLeaf)

			) %>%
			arrange(isAncestorLeaf) %>%
			mutate(my_C = 1:n()) %>%
			select(name, my_C)
	) %>%
		pull(my_C)
}

# @description Parse the stan fit object
parse_summary = function(fit) {
	fit %>%
		rstan::summary() %$% summary %>%
		as_tibble(rownames = ".variable") %>%
		filter(grepl("prop", .variable)) %>%
		separate(.variable,
						 c(".variable", "Q", "C"),
						 sep = "[\\[,\\]]",
						 extra = "drop") %>%
		mutate(C = C %>% as.integer, Q = Q %>% as.integer)
}

median_qi_nest_draws = function(d){
	# Anonymous function to add the draws to the summary

		left_join(
			d %>%
				#group_by(.variable,  Q,  C) %>%
				tidybayes::median_qi() %>%
				ungroup(),
			# Reattach draws as nested
			d %>%
				ungroup() %>%
				nest(.draws = c(.chain, .iteration, .draw , .value, .value_relative)),
			by = c(".variable", "Q", "sample", "C", "Cell type category", "level")
		)
}

# @description Parse the stan fit object and check for divergences
parse_summary_check_divergence = function(draws) {
	draws %>%

		group_by(level, .variable,  Q, sample,  C, `Cell type category`) %>%

		# If not converged choose the majority chains
		mutate(converged = diptest::dip.test(`.value_relative`) %$%	`p.value` > 0.05) %>%

		# If some proportions have not converged chose the most populated one
		# do(
		# 	(.) %>%
		# 		ifelse_pipe(
		# 			(.) %>% distinct(converged) %>% pull(1) %>% `!`,
		# 			~ .x %>% choose_chains_majority_roule
		# 		)
		# ) %>%

		# Anonymous function - add summary fit to converged label
		# input: tibble
		# output: tibble
		{
			left_join(
				(.) %>% select(-converged) %>% median_qi_nest_draws(),
				(.) %>% distinct(converged),
				by = c("level", ".variable", "Q", "sample", "C", "Cell type category")
			)
		} %>%
		ungroup()
}

# @description create tree object that is in data directory
create_tree_object = function(my_ref = ARMET::ARMET_ref) {
	

	
	#yaml:: yaml.load_file("~/PhD/deconvolution/ARMET/data/tree.yaml") %>%
	tree = 
		yaml::yaml.load_file("data/tree.yaml") %>%
		data.tree::as.Node() %>%

		{

			my_ct = my_ref %>% distinct(`Cell type category`) %>% pull(1) %>% as.character

			# Filter if not in referenc
			data.tree::Prune(., pruneFun = function(x) ( x$name %in% my_ct ))

			# Sort tree by name
			data.tree::Sort(., "name")

			# Add C indexes
			.$Set(C =
							tibble(
								name = .$Get('name'),
								level = .$Get('level')
							) %>%
							left_join((.) %>% arrange(level, name) %>%	 	mutate(C = 0:(n(

							) - 1)))	%>%
							pull(C))
			.$Set(C1 = get_idx_level(., 1))
			.$Set(C2 = get_idx_level(., 2))
			.$Set(C3 = get_idx_level(., 3))
			.$Set(C4 = get_idx_level(., 4))
			#		if(max(levels)>1) for(l in 2:max(levels)) { my_c = sprintf("C%s", l); .$Set(	my_c = get_idx_level(.,2)	); . }

			# Set Cell type category label
			.$Set("Cell type category" = .$Get("name"))

			.

		}
	
	save(tree, file="data/tree.rda", compress = "gzip")
	
	ancestor_child = tree %>% get_ancestor_child
	
	save(ancestor_child, file="data/ancestor_child.rda", compress = "gzip")
}

# @description Convert tibble to matrix
as_matrix = function(tbl, rownames = NULL) {
	# If matriix empty ski the whole thing
	if (length(tbl) == c(0))
		return(matrix()[0, 0])


	tbl %>%

		# Check if data frame is not numerical beside the rownames column (if present)
		{
			if (!tbl %>%
					{
						if (!is.null(rownames))
							(.) %>% dplyr::select(-contains(rownames))
						else
							(.)
					} %>%
					dplyr::summarise_all(class) %>%
					tidyr::gather(variable, my_class) %>%
					pull(my_class) %>% unique %>% identical("numeric"))
				warning("to_matrix says: there are NON-numerical columns, the matrix will NOT be numerical")

			(.)

		} %>%
		as.data.frame() %>%

		# Deal with rownames column if present
		{
			if (!is.null(rownames))
				(.) %>%
				magrittr::set_rownames(tbl %>% pull(!!rownames)) %>%
				dplyr::select(-!!rownames)
			else
				(.)
		} %>%

		# Convert to matrix
		as.matrix()
}

# @description Extension of data.tree package. It converts the tree into data frame

ToDataFrameTypeColFull = function(tree, fill = T, ...) {
	t = tree %>% data.tree::Clone()
	
	1:(t %$% Get("level") %>% max) %>%
		map_dfr(
			~ data.tree::Clone(t) %>%
				{
					data.tree::Prune(., function(x)
						x$level <= .x + 1)
					.
				} %>%
				data.tree::ToDataFrameTypeCol() %>%
				as_tibble
			
		) %>%
		distinct() %>%
		
		when(
			fill & ("level_3" %in% colnames(.)) ~ mutate(., level_3 = ifelse(level_3 %>% is.na, level_2, level_3)),
			TRUE ~ (.)
		) %>%
		when(
			fill & ("level_4" %in% colnames(.)) ~ mutate(., level_4 = ifelse(level_4 %>% is.na, level_3, level_4)),
			TRUE ~ (.)
		) %>%
		when(
			fill & ("level_5" %in% colnames(.)) ~ mutate(., level_5 = ifelse(level_5 %>% is.na, level_4, level_5)),
			TRUE ~ (.)
		) %>%
		when(
			fill & ("level_6" %in% colnames(.)) ~ mutate(., level_6 = ifelse(level_6 %>% is.na, level_5, level_6)),
			TRUE ~ (.)
		) %>%

		dplyr::select(..., everything())
	
}

#' vb_iterative
#' @keywords internal
#'
#' @description Runs iteratively variational bayes until it suceeds
#'
#' @importFrom rstan vb
#'
#' @param model A Stan model
#' @param output_samples An integer of how many samples from posteriors
#' @param iter An integer of how many max iterations
#' @param tol_rel_obj A real
#' @param algorithm A character
#'
#' @return A Stan fit object
#'
vb_iterative = function(model,
												output_samples,
												iter,
												tol_rel_obj,
												algorithm = "fullrank",
												...) {
	res = NULL
	i = 0
	while (res %>% is.null | i > 5) {
		res = tryCatch({
			my_res = vb(
				model,
				output_samples = output_samples,
				iter = iter,
				tol_rel_obj = tol_rel_obj,
				algorithm = algorithm,
				...
				#, pars=c("counts_rng", "exposure_rate", additional_parameters_to_save)
			)
			boolFalse <- T
			return(my_res)
		},
		error = function(e) {
			i = i + 1
			writeLines(sprintf("Further attempt with Variational Bayes: %s", e))
			return(NULL)
		},
		finally = {

		})
	}

	return(res)
}

get_level_lpdf_weights = function(df) {
	df %>%
		distinct(level, `Cell type category`) %>%
		drop_na %>%
		count(level) %>%
		mutate(`max n` = (.) %>% tail(1) %>% pull(n)) %>%
		mutate(weight = (`max n`) / n) %>%
		select(level, weight)
}

get_NB_qq_values = function(input.df, transcript_column) {
	transcript_column = enquo(transcript_column)

	input.df %>%
		group_by(!!transcript_column) %>%
		do({
			(.) %>%
				mutate(
					predicted_NB =

						qnbinom(
							# If 1 sample, just use median
							switch(
								((.) %>% nrow > 1) %>% `!` %>% `+` (1),
								ppoints(`read count normalised bayes`),
								0.5
							),
							size = .$sigma_inv_log %>% unique %>% exp %>% `^` (-1),
							mu = .$lambda_log %>% unique %>% exp
						)

				)
		}) %>%
		ungroup
}

filter_house_keeping_query_if_fixed =  function(.data, full_bayesian) {
	.data %>%

		# If full Bayesian false just keep house keeping
		ifelse_pipe(!full_bayesian,
								~ .x %>% filter(`house keeping` & `query`))
}


ref_mix_format = function(ref, mix) {
	bind_rows( 
		# Get reference based on mix genes
		ref %>% mutate(`query` = FALSE),
		
		# Select only markers
		mix %>%
			inner_join(ref %>% distinct(symbol, `house keeping`), by = "symbol") %>%
			mutate(`Cell type category` = "query") %>%
			mutate(`query` = TRUE)
	)	%>%

		# Add marker symbol indexes
		left_join((.) %>%
								filter(!`house keeping`) %>%
								distinct(symbol) %>%
								mutate(M = 1:n()), by = "symbol") %>%

		# Add sample indeces
		arrange(!`query`) %>% # query first
		mutate(S = factor(sample, levels = .$sample %>% unique) %>% as.integer) %>%

		# Add query samples indeces
		left_join((.) %>%
								filter(`query`) %>%
								distinct(sample) %>%
								mutate(Q = 1:n()), by="sample") %>%

		# Add house keeping into Cell type label
		mutate(`Cell type category` = ifelse(`house keeping`, "house_keeping", `Cell type category`)) %>%

		# Still needed?
		anti_join(
			(.) %>%
				filter(`house keeping` & !`query`) %>%
				distinct(symbol, level) %>%
				group_by(symbol) %>%
				arrange(level) %>%
				slice(2:max(n(), 2)) %>% # take away house keeping from level 2 above
				ungroup(),
			by = c("level", "symbol")
		) %>%

		# If house keeping delete level infomation
		mutate(level = ifelse(`house keeping`, NA, level)) %>%

		# Create unique symbol ID
		unite(ct_symbol, c("Cell type category", "symbol"), remove = F) %>%

		# Add gene idx
		left_join(
			(.) %>%
				filter(!`query`) %>%
				distinct(`Cell type category`, ct_symbol, `house keeping`) %>%
				arrange(!`house keeping`, ct_symbol) %>% # House keeping first
				mutate(G = 1:n()),
			by = c("ct_symbol", "Cell type category", "house keeping")
		) %>%
		left_join(
			(.) %>%
				filter(!`house keeping` & !`query`) %>%
				distinct(level, symbol) %>%
				arrange(level, symbol) %>%
				mutate(GM = 1:n()) %>%
				select(-level),
			by = "symbol"
		)

}

get_prop = function(fit, approximate_posterior, df, tree) {
	fit %>%

		# If MCMC is used check divergencies as well
		ifelse_pipe(
			!approximate_posterior,
			~ .x %>% parse_summary_check_divergence(),
			~ .x %>% parse_summary() %>% rename(.value = mean)
		) %>%

		# Parse
		separate(.variable, c(".variable", "level"), convert = T) %>%

		# Add tree information
		left_join(
			tree %>% data.tree::ToDataFrameTree("name", "C1", "C2", "C3", "C4") %>%
				as_tibble %>%
				select(-1) %>%
				rename(`Cell type category` = name) %>%
				gather(level, C,-`Cell type category`) %>%
				mutate(level = gsub("C", "", level)) %>%
				drop_na %>%
				mutate(C = C %>% as.integer, level = level %>% as.integer)
		) %>%

		# Add sample information
		left_join(df %>%
								filter(`query`) %>%
								distinct(Q, sample))
}

draws_to_alphas = function(.data, pars) {
	.data %>%
		tidybayes::gather_draws(`prop_[1,a-z]`[Q, C], regex = T) %>%
		ungroup() %>%
		# {
		# 	print((.))
		# 	(.)
		# } %>%
		filter(.variable %in% pars) %>%
		# {
		# 	print((.))
		# 	(.)
		# } %>%
		nest(data = -c(.variable, Q)) %>%
		mutate(
			alphas = map(
				data,
				~
					.x %>%
				  drop_na() %>%
				  spread(C, .value) %>%
					select(-c(1:3)) %>%
					as_matrix() %>%
					sirt::dirichlet.mle() %$%
					alpha %>%
					as_tibble() %>%
					rename(alpha = value) %>%
					mutate(C = 1:n()) %>%
				  spread(C, alpha)
			)
		) %>%
		select(-data, -Q) %>%
    unnest(alphas) %>%
    nest(alphas = -.variable) %>%
    mutate(alphas = map(alphas, ~ .x %>% select_if(function(x){!all(is.na(x))})))%>%
		pull(alphas)
}

get_null_prop_posterior = function(ct_in_nodes) {
	prop_posterior = list()
	for (i in 1:(length(ct_in_nodes))) {
		prop_posterior[[i]] =  matrix(ncol = ct_in_nodes[i]) %>% as_tibble() %>% setNames(c(1:ct_in_nodes[i])) %>% slice(0)
	}
	names(prop_posterior) = sprintf("prop_%s_prior", c(1, letters[1:(length(prop_posterior)-1)]))
	prop_posterior
}

plot_boxplot = function(input.df, symbols) {
	input.df %>%
		filter(symbol %in% symbols) %>%
		#filter(level == 1) %>%
		#filter(`Cell type category` == "endothelial")  %>%
		ggplot(
			aes(
				x = `Cell type category`,
				y = `read count normalised bayes` + 1,
				label = symbol,
				color = `regression`
			)
		) +
		geom_jitter() +
		geom_boxplot() +
		facet_wrap(~ symbol + `Cell type category`, scales = "free")  +
		expand_limits(y = 1, x = 1) +
		scale_y_log10()
}

plot_signatures = function(.data, n_markers, mix, ct1, ct2, n = 10, level) {
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
			axis.text.x = element_text(
				angle = 90,
				hjust = 1,
				vjust = 0.5
			),
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
			))
		)

	ARMET::ARMET_ref %>%
		left_join(n_markers, by = c("ct1", "ct2")) %>%
		filter_reference(mix) %>%
		filter(level %in% level) %>%
		filter(ct1 == !!ct1 & ct2 == !!ct2) %>%
		filter(rank < n) %>%
		distinct(
			`Cell type category`,
			sample,
			`read count normalised bayes`,
			symbol,
			regression,
			bimodality_NB
		) %>%
		ggplot(
			aes(
				x = `Cell type category`,
				y = `read count normalised bayes` + 1,
				label = symbol,
				color = `bimodality_NB`
			)
		) +
		geom_boxplot() +
		geom_jitter() +
		facet_wrap(~ symbol , scales = "free")  +
		expand_limits(y = 1, x = 1) +
		scale_y_log10() +
		my_theme
}

#' This is a generalisation of ifelse that acceots an object and return an objects
#' @keywords internal
#'
#' @import dplyr
#' @importFrom purrr as_mapper
#'
#' @param .x A tibble
#' @param .p A boolean
#' @param .f1 A function
#' @param .f2 A function
#'
#'
#' @return A tibble
ifelse_pipe = function(.x, .p, .f1, .f2 = NULL) {
	switch(.p %>% `!` %>% sum(1),
				 as_mapper(.f1)(.x),
				 if (.f2 %>% is.null %>% `!`)
				 	as_mapper(.f2)(.x)
				 else
				 	.x)

}

# This is a generalisation of ifelse that acceots an object and return an objects
level_to_plot_inferred_vs_observed  = function(result, level, S = NULL, cores = 20){

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
			axis.text.x = element_text(
				angle = 90,
				hjust = 1,
				vjust = 0.5
			),
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
			))
		)

	result$signatures[[level]] %>%
		filter(!query & !`house keeping`) %>%
		distinct(`Cell type category`, C, level, G, GM, symbol, lambda_log, sigma_inv_log) %>%
		{ print(1); Sys.time(); (.) } %>%

		# Add divergence
		left_join(
			result$fit[[level]] %>% rstan::summary() %$% summary %>% as_tibble(rownames = "par") %>% filter(Rhat > 1.6) %>%
				filter(grepl("^lambda_log", par)) %>%
				separate(par, c("par", "G"), sep="\\[|\\]", extra = "drop") %>%
				distinct(G) %>% mutate(G = G %>% as.integer) %>%
				mutate(converged = F)
		) %>%
		mutate(converged = ifelse(converged %>% is.na, T, converged)) %>%
		{ print(2); Sys.time(); (.) } %>%

		# If inferred replace lambda_log and sigma_inv_log
		ifelse_pipe(
			result$input$full_bayesian,
			~ .x %>% select(-lambda_log, -sigma_inv_log) %>%
				left_join(
					result$fit[[level]] %>% tidybayes::spread_draws(lambda_log[G], sigma_inv_log[G]) %>% ungroup() %>%
						rename(lambda_log = lambda_log,sigma_inv_log = sigma_inv_log)
				)
		) %>%
		{ print(3); Sys.time(); (.) } %>%

		left_join(
			ARMET::ARMET_ref %>% distinct(symbol, ct1, ct2, level)
		) %>%

		# Add proportions
		left_join(
			result$proportions %>% select(level, Q, sample, .draws, `Cell type category`) %>% unnest(.draws)
		) %>%

		# add expsure
		left_join(
			result$fit[[level]] %>% tidybayes::spread_draws(exposure_rate[S]) %>% ungroup() %>% rename(Q = S)
		) %>%

		# Filter by sample
		ifelse_pipe(
			S %>% is.null %>% `!`,
			~ .x %>% filter(Q == S)
		)	 %>%

		# Calculate sum
		mutate(
			lambda_exp = lambda_log %>% exp,
			sigma_exp = 1 / exp(sigma_inv_log)
		) %>%

		{ print(4); Sys.time(); (.) } %>%

		# Filter just first 30 draws
		inner_join( (.) %>% distinct(.draw) %>% sample_n(30) ) %>%

		do_parallel_start(cores, "symbol") %>%
		do({

			`%>%` = magrittr::`%>%`

			sum_NB = function(lambda, sigma, prop){

				prop_mat = matrix(prop, nrow=1)
				lambda_mat = matrix(lambda, ncol = 1)
				sigma_mat = matrix(sigma, ncol = 1)

				lambda_sum = prop_mat %*% lambda_mat;
				sigma_sum =
					lambda_sum^2 /
					(
						prop_mat %*%
							(
								lambda_mat^2 /
									sigma_mat
							)
					) ;


				c(lambda_sum, sigma_sum)

			}

			(.) %>%
				tidyr::nest(data_for_sum = -c(level, symbol, GM, converged, ct1    ,     ct2   ,     sample   ,   exposure_rate,   .chain, .iteration ,.draw )) %>%
				dplyr::mutate(.sum = purrr::map(
					data_for_sum,
					~
						sum_NB(.x$lambda_exp, .x$sigma_exp, .x$.value) %>%
						as.matrix %>%
						t %>%
						tibble::as_tibble() %>%
						setNames(c("lambda_sum", "sigma_sum"))

				))

		}) %>%
		do_parallel_end() %>%

		{ print(5); Sys.time(); (.) } %>%

		select(-data_for_sum) %>%
		unnest(.sum) %>%
		{ print(6); Sys.time(); (.) } %>%

		# normalise
		mutate(lambda_sum = lambda_sum * exp(exposure_rate)) %>%

		# Calculate generated quantities
		mutate(counts_inferred = rnbinom(n(), mu = lambda_sum, size = sigma_sum)) %>%
		{ print(7); Sys.time(); (.) } %>%

		# Summarise
		group_by(level, symbol, GM, converged, ct1, ct2, sample, .chain) %>%
		summarize(.mean = mean(counts_inferred),
							.sd = sd(counts_inferred),
							.q025 = quantile(counts_inferred, probs = .025),
							.q25 = quantile(counts_inferred, probs = .25),
							.q50 = quantile(counts_inferred, probs = .5),
							.q75 = quantile(counts_inferred, probs = .75),
							.q97.5 = quantile(counts_inferred, probs = .975)
		) %>%
		{ print(8); Sys.time(); (.) } %>%

		# Add counts
		left_join(	result$input$mix %>%	gather(symbol, count, -sample) ) %>%
		{ print(9); Sys.time(); (.) } %>%

	 { ((.) %>%	ggplot(aes(x = count+1, y=.q50+1, color=ct1, shape = converged, GM = GM)) +
		geom_abline() +
		geom_errorbar(aes(ymin = .q025 + 1, ymax=.q97.5 + 1), alpha=0.5) +
		geom_point() +
		scale_y_log10() +
		scale_x_log10() +
		facet_grid(converged ~.chain) +
		my_theme) %>% print

	 	(.)
	 }
}

# level_to_plot_inferred_vs_observed(result, 3)

#' Create the design matrix
#' @keywords internal
#'
#' @param input.df A tibble
#' @param formula A formula
#' @param sample_column A symbol
#' 
create_design_matrix = function(input.df, formula, sample_column){

	sample_column = enquo(sample_column)

	model.matrix(
		object = formula,
		data =
			input.df %>%
			select(!!sample_column, one_of(parse_formula(formula)$covariates)) %>%
			distinct %>% arrange(!!sample_column)

	)

}

# Formula parser
parse_formula <- function(fm) {
	
	components = as.character(attr(terms(fm), "variables"))[-1]

	components_formatted = 
		components %>% 
		map_chr(~ .x %>% 
							gsub("censored\\(|\\)| ", "", .) %>% 
							str_split("\\,") %>% .[[1]] %>% .[1]
		) 
	
	covariates_formatted = components_formatted %>% when(attr(terms(fm), "response") == 1 ~ (.)[-1], (.))
	response_formatted = components_formatted %>% when(attr(terms(fm), "response") == 1 ~ (.)[1])
	
	censored_formatted= 
		components %>% 
		grep("censored", ., value = T) %>% 
		map_chr(~ .x %>% 
							gsub("censored\\(|\\)| ", "", .) %>% 
							str_split("\\,") %>% .[[1]] %>% .[1]
					) 
	
	censored_column = 
		components %>% 
		grep("censored(", ., fixed = T, value = T)  %>% 
		gsub("censored\\(|\\)| ", "", .) %>% str_split("\\,") %>% 
		when(length(.)>0 ~(.) %>% .[[1]] %>% .[-1])
	
	censored_value_column = 
		components%>% grep("censored(", ., fixed = T, value = T)  %>% 
		gsub("censored\\(|\\)| ", "", .) %>% str_split("\\,") %>% 
		when(length(.)>0 ~(.) %>% .[[1]]%>% .[1])
	
	formula_formatted =
		fm %>%
		when(
			length(grep( "censored", components, value = T)) == 0 ~ fm,
			~ Reduce(paste, deparse(fm)) %>%
				gsub( grep( "censored", components, value = T), censored_formatted, ., fixed = T)	 %>%
				as.formula
		)
		
	
	formula_censored_formatted =
		fm %>%
		when(
			length(grep( "censored", components, value = T)) == 0 ~ NULL,
			~ Reduce(paste, deparse(fm)) %>%
				gsub( grep( "censored", components, value = T), "proportion", ., fixed = T)	 %>%
				sprintf("%s%s", censored_formatted, .) %>%
				as.formula
		)
	
	
	list(
		components = components,
		components_formatted = components_formatted,
		covariates_formatted = covariates_formatted,
		response_formatted = response_formatted,
		censored_formatted = censored_formatted,
		censored_column = censored_column,
		censored_value_column = censored_value_column,
		formula_formatted = formula_formatted,
		formula_censored_formatted = formula_censored_formatted
	)
}

rebuild_last_component_sum_to_zero = function(.){
	 
	(.) %>%
		nest(data = -c(.variable, A)) %>%
		mutate(data = map(data, ~.x %>%
												mutate(C = C +1) %>%
												
												# Add a 0 component, the first one
												bind_rows({
													
													(.) %>%
														
														# Select one of the C does not matter
														filter(C ==2) %>%
														mutate(C = rep(1, n())) %>%
														mutate(.value = rep(0, n())) 
													
												})  %>%
												group_by(.chain ,.iteration, .draw) %>%
												mutate(.value = .value - mean(.value)) %>%
												ungroup  %>%
												arrange(C)
											
		)) %>%
		unnest(data)
	
}

get_relative_zero = function(fit_parsed){

	fit_parsed %>%
		nest(data = -.variable) %>%
		mutate(zero = map(
			data,
			~ {

				rng = .x %>% pull(.value) %>% summary %>% `[` (c(1,6))

				.x %>%
				filter(A == 2) %>%
				nest(data2 = -C) %>%
				mutate(zero_C = map(data2, ~ .x %>%
															pull(.value) %>%
															density(na.rm = T, from = rng[1], to =rng[2]) %>%
															{	tibble(x = (.)$x, y = (.)$y)	} %>%
															mutate(y = y/max(y))
															)) %>%
				unnest(zero_C) %>%
				select(-data2) %>%
				group_by(x) %>%
				summarise(y = sum(y)) %>%
				arrange(y %>% desc) %>% slice(1) %>% select(zero = x)
		})) %>%
		unnest(cols = c(data, zero))

}

#' Get column names either from user or from attributes
#' @keywords internal
#'
#' @importFrom rlang quo_is_symbol
#'
#' @param .data A tibble
#' @param .sample A character name of the sample column
#' @param .transcript A character name of the transcript/gene column
#' @param .abundance A character name of the read count column
#'
#' @return A list of column enquo or error
get_sample_transcript_counts = function(.data, .sample, .transcript, .abundance){


	my_stop = function() {
		stop("
        ARMET says: The fucntion does not know what your sample, transcript and counts columns are.\n
        You have to either enter those as symbols (e.g., `sample`), \n
        or use the funtion create_tt_from_tibble() to pass your column names that will be remembered.
      ")
	}

	if( .sample %>% quo_is_symbol() ) .sample = .sample
	else my_stop()

	if( .transcript %>% quo_is_symbol() ) .transcript = .transcript
	else my_stop()

	if( .abundance %>% quo_is_symbol() ) .abundance = .abundance
	else my_stop()

	list(.sample = .sample, .transcript = .transcript, .abundance = .abundance)

}

eliminate_sparse_transcripts = function(.data, .transcript){
	# Parse column names
	.transcript = enquo(.transcript)

	warning("Some transcripts have been omitted from the analysis because not present in every sample.")

	.data %>%
		add_count(symbol, name = "my_n") %>%
		filter(my_n == max(my_n)) %>%
		select(-my_n)
}

get_specific_annotation_columns = function(.data, .col){
	
	
	# Comply with CRAN NOTES
	. = NULL
	
	# Make col names
	.col = enquo(.col)
	
	# x-annotation df
	n_x = .data %>% distinct(!!.col) %>% nrow
	
	# Sample wise columns
	.data %>%
		select(-!!.col) %>%
		colnames %>%
		map(
			~
				.x %>%
				ifelse_pipe(
					.data %>%
						distinct(!!.col, !!as.symbol(.x)) %>%
						nrow %>%
						equals(n_x),
					~ .x,
					~ NULL
				)
		) %>%
		
		# Drop NULL
		{	(.)[lengths((.)) != 0]	} %>%
		unlist
	
}

#' @importFrom tidybulk as_matrix
#' @keywords internal
#' 
#' @param fit_parsed A fit object
prop_to_list = function(fit_parsed){
	fit_parsed %>%
		median_qi() %>%
		drop_na  %>%
		ungroup() %>% 
		nest(data = -.variable) %>% 
		mutate(data =
					 	map(
					 		data, 
					 		~.x %>% 
					 			select(Q, C, .value) %>% 
					 			spread(C, .value) %>% 
					 			tidybulk::as_matrix(rownames = "Q")
					 	)
					) %>%
		{ x = (.); x %>% pull(data) %>% setNames(x %>% pull( .variable))}
}

gamma_alpha_beta = function(x){
	summ = summary((dglm(x~1, family=Gamma(link="log"), mustart=mean(x))))
	mu <- exp(summ$coefficients[1])
	shape <- exp(-summ$dispersion)
	scale <- mu/shape
	c(shape, 1/scale)
}

get_ancestor_child = function(tree){
	tree %>% ToDataFrameTypeColFull %>% distinct(level_1, level_2) %>% setNames(c("ancestor", "Cell type category")) %>% bind_rows(
		tree %>% ToDataFrameTypeColFull %>% distinct(level_2, level_3) %>% setNames(c("ancestor", "Cell type category"))
	) %>%
		bind_rows(
			tree %>% ToDataFrameTypeColFull %>% distinct(level_3, level_4) %>% setNames(c("ancestor", "Cell type category"))
		) %>%
		bind_rows(
			tree %>% ToDataFrameTypeColFull %>% distinct(level_4, level_5) %>% setNames(c("ancestor", "Cell type category"))
		) %>%
		filter(ancestor != `Cell type category`)
}

#' @importFrom purrr map_int
#' @importFrom foreach foreach
#' @importFrom foreach %do%
#' @keywords internal
#' 
#' @param tree A tree
get_tree_properties = function(tree){
	

	# Set up tree structure
	levels_in_the_tree = 1:4
	
	ct_in_nodes =
		tree %>%
		data.tree::ToDataFrameTree("name", "level", "C", "count", "isLeaf") %>%
		as_tibble %>%
		arrange(level, C) %>%
		filter(!isLeaf) %>%
		pull(count)
	
	ct_in_levels = foreach(l = levels_in_the_tree + 1, .combine = c) %do% {
		data.tree::Clone(tree) %>%
			ifelse_pipe((.) %>% data.tree::ToDataFrameTree("level") %>% pull(2) %>% max %>% `>` (l),
									~ {
										.x
										data.tree::Prune(.x, function(x)
											x$level <= l)
										.x
									})  %>%
			data.tree::Traverse(., filterFun = isLeaf) %>%
			length()
	}
	
	# # Get the number of leafs for every level
	# ct_in_levels =
	# 	levels_in_the_tree + 1 %>%
	# 	map_int(~ {
	# 		data.tree::Clone(tree) %>%
	# 			when(
	# 				data.tree::ToDataFrameTree(., "level") %>% pull(2) %>% max %>% gt(.x) ~ {
	# 					t = (.)
	# 					data.tree::Prune(t, function(x)	x$level <= .x)
	# 					t
	# 				},
	# 				~ (.)
	# 			) %>%
	# 			data.tree::Traverse(., filterFun = isLeaf) %>%
	# 			length()
	# 	})
	
	n_nodes = ct_in_nodes %>% length
	n_levels = ct_in_levels %>% length
	
	# Needed in the model
	singles_lv2 = tree$Get("C1", filterFun = isLeaf) %>% na.omit %>% as.array
	SLV2 = length(singles_lv2)
	parents_lv2 = tree$Get("C1", filterFun = isNotLeaf) %>% na.omit %>% as.array
	PLV2 = length(parents_lv2)
	
	singles_lv3 = tree$Get("C2", filterFun = isLeaf) %>% na.omit %>% as.array
	SLV3 = length(singles_lv3)
	parents_lv3 = tree$Get("C2", filterFun = isNotLeaf) %>% na.omit %>% as.array
	PLV3 = length(parents_lv3)
	
	singles_lv4 = tree$Get("C3", filterFun = isLeaf) %>% na.omit %>% as.array
	SLV4 = length(singles_lv4)
	parents_lv4 = tree$Get("C3", filterFun = isNotLeaf) %>% na.omit %>% as.array
	PLV4 = length(parents_lv4)
	
	list(
		ct_in_nodes =ct_in_nodes,
		
		# Get the number of leafs for every level
		ct_in_levels = ct_in_levels,
		
		n_nodes = ct_in_nodes %>% length,
		n_levels = ct_in_levels %>% length,
		
		# Needed in the model
		singles_lv2 = singles_lv2,
		SLV2 = SLV2,
		parents_lv2 = parents_lv2,
		PLV2 = PLV2,
		
		singles_lv3 = singles_lv3,
		SLV3 = SLV3,
		parents_lv3 = parents_lv3,
		PLV3 = PLV3,
		
		singles_lv4 = singles_lv4,
		SLV4 = SLV4,
		parents_lv4 = parents_lv4,
		PLV4 = PLV4
	)
}

#' @importFrom tidygraph tbl_graph
#' @keywords internal
#' 
#' @param .data A tibble
#' @param credible_interval A double
cluster_posterior_slopes = function(.data, credible_interval = 0.67){
	   
	# Add cluster info to cell types per node
	.data %>%
		filter(.variable %>% is.na %>% `!`) %>%
		nest(node = -c(level, .variable)) %>%
		mutate(node = map(
			node,
			~ {
				.x %>% 
					left_join(
						
						# Unnest data
						(.) %>% 
							extract_CI(credible_interval) %>%
							select(-c(proportions  ,   draws , rng_prop )) %>%
							
							# Build combination of cell types
							combine_nest(
								.names_from = `Cell type category`,
								.values_from = c(.lower_alpha2, .upper_alpha2, Rhat)
							)  %>%
							
							# Check overlap of credible intervals
							mutate(is_cluster = map_lgl(
								data,
								~ .x %>% 
									summarise(ma = max(.lower_alpha2), mi = min(.upper_alpha2), converged = any(Rhat > 1.6 & !is.na(Rhat))==F) %>% 
									mutate(is_cluster = ma < mi & converged) %>%
									pull(is_cluster)
							)) %>% 
							
							# If there is not cluster
							
							# Find communities based on cell type clusters
							{
								ct_levels = (.) %>% arrange(!is_cluster) %>% select(1:2) %>% as_matrix %>% t %>% as.character() %>% unique
								
								(.) %>%
									filter(is_cluster) %>% 
									select(1:2) %>%
										tbl_graph(
										edges = .,
										nodes = data.frame(name = ct_levels)
									)%>%
									mutate(community = as.factor(tidygraph::group_infomap())) 
							}	%>%
							
							# Format for joining
							as_tibble() %>%
							rename(`Cell type category` = name),
						by = "Cell type category"
					)
				
			}	)) %>%
		unnest(node)
	
}

identify_baseline_by_clustering = function(.data, CI){
	  
	.data %>%
		nest(node = -.variable) %>%
		
		mutate(node = map(
			node, ~ .x %>%
				#mutate(baseline = TRUE) %>%
				add_count(community) %>%
				
				# Add sd
				nest(comm_data = -community) %>%
				mutate(
					sd_community = map_dbl( comm_data, ~ .x %>% unnest(draws) %>% filter(A ==2) %>% pull(.value) %>% sd ),
					median_community = map_dbl( comm_data, ~ .x %>% unnest(draws) %>% filter(A ==2) %>% pull(.value) %>% median )
				) %>%
				unnest(comm_data) %>%
				
				purrr::when(
					
					# # If I have a reference change make zero frm it
					# (.) %>% distinct(fold_change_ancestor) %>% pull(1) %>% `!=` (0) ~ 
					# 	(.) %>% mutate(zero = -fold_change_ancestor),
					
					# Otherwise, if I have a consensus of overlapping posteriors make that zero
					(.) %>% distinct(community, n) %>% count(n) %>% nrow %>% `>` (1) ~ 
						(.) %>% mutate(zero = (.) %>% filter( n == max(n)) %>% pull(median_community) %>% mean),
					
					# Otherwise make zero the 0
					~ (.) %>% mutate(zero = 0)
				)
		)) %>%
				 
		# 		
		# 		# If we have a unique bigger community
		# 		ifelse2_pipe(
		# 			(.) %>% distinct(community) %>% nrow %>% equals(1),
		# 			(.) %>% distinct(community, n) %>% count(n) %>% arrange(n %>% desc) %>% slice(1) %>% pull(nn) %>% equals(1) ,
		# 			
		# 			# If I have just one community
		# 			~ .x %>% mutate(baseline = TRUE),
		# 			
		# 			# Majority roule
		# 			~ .x %>% mutate(baseline = n == max(n)) ,
		# 			
		# 			# If ancestor changed, if no consensus the zero will be absolute 0 
		# 			~ .x %>% mutate(baseline = FALSE)
		# 		)
		# 	
		# )) %>% 
		# 
		# # Select zero. If I hav comunity select mean otherwise select 0
		# mutate(zero = map_dbl(node, ~ .x %>% filter(baseline) %>% ifelse_pipe( (.) %>% distinct(community) %>% nrow %>% equals(1), ~.x %>% pull(mean_community) %>% unique, ~ 0) )) %>%
		unnest(node)
	
	
}

extract_CI =  function(.data, credible_interval = 0.90){
	 
		.data %>%
		mutate(regression = map(draws,
														~ .x %>%
															group_by(A) %>%
															tidybayes::median_qi(.width = credible_interval) %>%
															ungroup()  %>%
															pivot_wider(
																names_from = A,
																values_from = c(.value, .lower, .upper),
																names_prefix = "alpha"
															) 
		)) %>%
		unnest(cols = c(regression)) %>%
		
		# If present the second test 
		when("draws_cens" %in% colnames(.) ~ (.) %>% left_join(
			.data %>%
				mutate(regression = map(draws_cens,
																~ .x %>%
																	group_by(A) %>% 
																	mutate(.draw = 1:n()) %>%
																	select(-one_of(".draw2")) %>%
																	
																	nest(draws = -c(A)) %>%
																	
																	mutate(prob_non_0 = map_dbl(draws, ~.x %>% draws_to_prob_non_zero)) %>%
																	mutate(draws = map(draws, ~.x %>% tidybayes::mean_qi())) %>%
																	unnest(draws) %>%
																	
																	
																	#group_by(.chain ,.iteration, .draw) %>%
																	# tidybayes::median_qi(.width = credible_interval) %>%
																	# ungroup() %>%
																	pivot_wider(names_from = A, values_from=c(.value, .value.lower, .value.upper,prob_non_0 ))
				)) %>%
				unnest(cols = c(regression)) %>%
				select(C, `Cell type category`,starts_with(".value"), ,starts_with(".value.upper"), ,starts_with(".value.lower"), starts_with("prob_non_0"))
			
		), 
		~ (.))
	
}

calculate_x_for_polar = function(.data){
	# Create annotation
	internal_branch_length = 40
	external_branch_length = 10
	
	
	# Integrate data
	tree_df_source = 
		ARMET::tree %>%
		data.tree::ToDataFrameTree("name", "isLeaf", "level", "leafCount", pruneFun = function(x)	x$level <= 5) %>%
		as_tibble() %>%
		rename(`Cell type category` = name) %>%
		mutate(level = level -1)
	
	# Calculate x
	map_df(
		0:4,
		~ tree_df_source %>%
			filter(level == .x | (level < .x & isLeaf)) %>%
			mutate(leafCount_norm = leafCount/sum(leafCount)) %>%
			mutate(leafCount_norm_cum = cumsum(leafCount_norm)) %>%
			mutate(length_error_bar = leafCount_norm - 0.005) %>%
			mutate(x = leafCount_norm_cum - 0.5 * leafCount_norm) 	
	) %>%
		
		# Attach data
		left_join(.data) %>%
		
		# process
		mutate(Estimate = ifelse(significant, fold_change, NA)) %>%
		mutate(branch_length = ifelse(isLeaf, 0.1, 2)) %>%
		
		# Correct branch length
		mutate(branch_length = ifelse(!isLeaf, internal_branch_length,	external_branch_length) ) %>%
		mutate(branch_length =  ifelse(	isLeaf, branch_length + ((max(level) - level) * internal_branch_length),	branch_length	)) 
	
}

get_draws = function(fit_prop_parsed, level, internals){
	
	my_c = as.symbol(sprintf("C%s", level))
	ancestor_c = as.symbol(sprintf("C%s", level-1))
	
	my_value_column = as.symbol(sprintf(".value%s", level) )
	ancestor_value_column = as.symbol(sprintf(".value%s", level-1) )

	internals$draws[[level-1]] %>%
		# Temporary I have to fix!!!
		when(level ==2 ~ (.) %>% rename(C1 = C, !!ancestor_value_column := .value), ~ (.)) %>%
		#
		left_join(
			fit_prop_parsed %>%
				left_join(
					######## ALTERED WITH TREE
					tibble(
						.variable = (.) %>% distinct(.variable) %>% pull(),
						!!ancestor_c := internals$tree_properties[[sprintf("parents_lv%s", level)]]
					)
					#
				) %>%
				select(-.variable) %>%
				rename(!!my_value_column := .value, !!my_c  := C) ,
			by = c( ".chain", ".iteration", ".draw", "Q", sprintf("C%s", level -1) )
		) %>%
		group_by(.chain, .iteration, .draw, Q) %>%
		arrange_at(vars(contains("C1"), contains("C2"), contains("C3"), contains("C4"), contains("C5"))) %>% 
		mutate(
			!!my_c := tree$Get(sprintf("C%s", level)) %>% na.omit,
			`Cell type category` = tree$Get(sprintf("C%s", level)) %>% na.omit %>% names
		) %>%
		ungroup() %>%
		mutate(.value_relative := !!my_value_column) %>%
		mutate(!!my_value_column := ifelse(!!my_value_column %>% is.na, !!ancestor_value_column, !!ancestor_value_column * !!my_value_column))
}

get_props = function(draws, level, df, approximate_posterior){
	
	my_c = as.symbol(sprintf("C%s", level))

	my_value_column = as.symbol(sprintf(".value%s", level) )

	draws %>%
		#when(level ==1 ~ (.) %>% rename(C1 = C), ~ (.)) %>%
		select(.chain,
					 .iteration,
					 .draw,
					 Q,
					 !!my_c ,
					 `Cell type category`,
					 !!my_value_column,
					 .value_relative) %>%
		rename(C := !!my_c, .value = !!my_value_column) %>%
		mutate(.variable = sprintf("prop_%s", level)) %>%
		mutate(level := !!level) %>%
		
		# add sample annotation
		left_join(df %>% distinct(Q, sample), by = "Q")	%>%
		
		# If MCMC is used check divergences as well
		ifelse_pipe(
			!approximate_posterior,
			~ .x %>% parse_summary_check_divergence(),
			~ .x %>% parse_summary() %>% rename(.value = mean)
		) %>%
		
		left_join(df %>% distinct(Q, sample))
	
}

get_alpha = function(fit, level){
	
	my_c = as.symbol(sprintf("C%s", level))
	
	fit %>%
		draws_to_tibble("alpha_", "A", "C") %>%

		# rebuild the last component sum-to-zero
		rebuild_last_component_sum_to_zero() %>%
		

		arrange(.chain, .iteration, .draw,     A) %>%
		
		nest(draws = -c(C, .variable)) %>%
		
		# Attach convergence information
		left_join(
			fit %>% summary_to_tibble("alpha_", "A", "C") %>% filter(A == 2) %>% 
				select(.variable, C, Rhat),
			by = c(".variable", "C")
		) %>%
		
		# FOR HIERARCHICAL
		mutate(C = 1:n()) %>%
		
		left_join(
			
			tree %>%
				ToDataFrameTree("name", "level", sprintf("C%s", level)) %>%
				filter(level == !!level+1) %>%
				arrange(!!my_c) %>%
				mutate(C = 1:n()) %>%
				select(name, C) %>%
				rename(`Cell type category` = name)
			
		) %>%
		
		# Attach generated quantities
		separate(.variable, c("par", "node"), remove = F) %>%
		# left_join(
		# 	fit %>% 
		# 		draws_to_tibble("prop_", "Q", "C") %>%
		# 		filter(grepl("_rng", .variable)) %>% 
		# 		mutate(Q = as.integer(Q)) %>%
		# 		mutate(.variable = gsub("_rng", "", .variable)) %>% 
		# 		separate(.variable, c("par", "node"), remove = F)  %>% 
		# 		select(-par) %>% drop_na() %>% nest(rng = -c(node, C)) %>%
		# 		mutate(C = 1:n()) 
		# ) %>%
		
		# Add level label
		mutate(level = !!level)

}

draws_to_tibble = function(fit, par, x, y) {
	 
	par_names = names(fit) %>% grep(sprintf("%s", par), ., value = T)
	
	fit %>%
		rstan::extract(par_names, permuted=F) %>% 
		as.data.frame %>% 
			as_tibble() %>%
			mutate(.iteration = 1:n()) %>% 
			pivot_longer(
				names_to = c("dummy", ".chain", ".variable", x, y),  
				cols = contains(par), 
				names_sep = "\\.|\\[|,|\\]|:",
				values_to = ".value"
			) %>%
		mutate(.chain = as.integer(.chain), !!as.symbol(x) := as.integer(!!as.symbol(x)), !!as.symbol(y) := as.integer(!!as.symbol(y))) %>%
		select(-dummy) %>%
			arrange(.variable, !!as.symbol(x), !!as.symbol(y), .chain) %>%
		group_by(.variable, !!as.symbol(x), !!as.symbol(y)) %>%
		mutate(.draw = 1:n()) %>%
		ungroup() %>%
		select(!!as.symbol(x), !!as.symbol(y), .chain, .iteration, .draw ,.variable ,     .value)
		
}

summary_to_tibble = function(fit, par, x, y) {
	
	par_names = names(fit) %>% grep(sprintf("%s", par), ., value = T)
	
	fit %>%
		rstan::summary(par_names) %$%
		summary %>%
		as_tibble(rownames = ".variable") %>% tidyr::extract(col = .variable, into = c(".variable", x, y), "(.+)\\[(.+),(.+)\\]", convert = T) 
	
	
}

get_alpha_test = function(slope, which_changing, cell_types){

	# Get the alpha matrix

	intercept = rep(0, length(cell_types))
	slope_arr = rep(0, length(cell_types))

	slope_arr[which_changing] = slope
	matrix(intercept %>%	c(slope_arr), ncol = 2)

}

get_survival_X = function(S){
	readRDS("dev/PFI_all_cancers.rds") %>%
		filter(PFI.2 == 1 & !is.na(PFI.time.2)) %>%
		select(real_days = PFI.time.2 ) %>%
		mutate(real_days = real_days %>% scale(center = F) %>% as.numeric) %>%
		sample_n(S) %>%
		mutate(sample = sprintf("S%s", 1:n())) %>%
		mutate(alive = sample(0:1, n(), replace = T)) %>%
		mutate(days = ifelse(alive==1, real_days/2, real_days) ) %>%
		mutate(intercept = 1)
}

#' @keywords internal
#' 
#' @importFrom nanny subset
#' @importFrom gtools rdirichlet 
#' 
#' @param .data A tibble
#' @param X_df A design matrix
#' @param alpha A real
#'
generate_mixture = function(.data, X_df, alpha) {
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
	
	
	X = X_df %>% select(intercept, real_days) %>% nanny::as_matrix()
	
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
		ct_changing %>% when(
			length(.) > 0 ~ 	matrix(c(rep(1, 2), c(0, 1)), ncol = 2)  %*% t(alpha) %>%
				apply(1, softmax) %>%
				t %>%
				`*` (40) %>%
				apply(1, softmax) %>%
				.[ct_names == ct_changing,] %>%
				{	max(.) / min(.)	} %>%
				{ slope = alpha[,2][ alpha[,2]!=0]; ifelse(slope<0, -(.), (.)) },
			~ 0
												 )
	
	# Add counts
	dirichlet_source =
		cell_type_proportions %>%
		left_join(.data, by = c("Cell type category", "sample"))
	
	# Make mix
	dirichlet_source %>%
		mutate(c = `count normalised bayes` * p) %>%
		group_by(run, symbol) %>%
		summarise(`count mix` = c %>% sum) %>%
		ungroup %>%
		
		left_join(dirichlet_source %>% nanny::subset(run) ) %>%
		
		mutate(fold_change = fold_change) %>%
		
		# Add proportions
		add_attr(cell_type_proportions, "proportions") 
	
}

get_noiseless_harmonised = function(){
	
	mix_base_unharmonized = readRDS("dev/mix_base_noiseless.RDS")
	
	my_markers =
		ARMET::ARMET_ref %>%
		
		left_join(ARMET::n_markers, by = c("ct1", "ct2")) %>%
		filter_reference(
			mix_base_unharmonized %>%
				filter(level == 3) %>%
				distinct(`Cell type category`, symbol, `count normalised bayes`) %>%
				spread(symbol, `count normalised bayes`),
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

run_censored_model_iterative_OLD = function(.data) {
	res = NULL
	i = 0
	while (res %>% is.null | i > 5) {
		res = tryCatch({
			my_res =	rstan::optimizing(
				stanmodels$censored_regression,
				data = .data) %$%
					
					# Formate results
					par %>%
					enframe() %>%
					filter(grepl("alpha", name)) %>%
					tidyr::extract(name, c("A", "C"), ".+\\[([0-9]),([0-9])\\]", convert = TRUE) %>%
					filter(A ==2) %>%
					mutate(C = !!.data[["colnames_prop_logit_scaled"]] %>% as.integer) 
			boolFalse <- T
			return(my_res)
		},
		error = function(e) {
			i = i + 1
			writeLines(sprintf("Further attempt with optimise: %s", e))
			return(NULL)
		},
		finally = {
		})
	}
	
	return(res)
}

run_censored_model_iterative = function(.data){
	
	sampling_iter = 5
	
	rstan::sampling(
		stanmodels$censored_regression,
		data = .data, chain = 1, iter=150+sampling_iter, warmup=150, save_warmup=F, refresh = 2000, init="0"
	)   %>%
		tidybayes::gather_draws(alpha[A, C]) %>%
		ungroup() %>%
		rename(value = .value, .draw2 = .draw) %>%
		select(-.variable, -.chain, -.iteration) %>%
	
		# Change C 
		nest(data = -C) %>%
		mutate(C = !!.data$prop_C_names %>% as.integer)  %>%
		unnest(data) 
	
}

run_censored_model = function(.data, sampling = F){
	 
	if(sampling)
		rstan::sampling(
			stanmodels$censored_regression,
			data = .data) %>%
		tidybayes::gather_draws(alpha[A, C]) %>%
	#	filter(A==dim(.data$X[1]))  %>%
		nest(draws = -c(A, C)) %>%
		
		mutate(prob_non_0 = map_dbl(draws, ~.x %>% draws_to_prob_non_zero)) %>%
		mutate(draws = map(draws, ~.x %>% tidybayes::mean_qi())) %>%
		unnest(draws) %>%
		
		#mutate(C = !!.data[["colnames_prop_logit_scaled"]] %>% as.integer)  %>%
		rename(value = .value)
	else
		run_censored_model_iterative(.data)
	
}

#' @keywords internal
#' @importFrom abind abind
#' @importFrom boot logit
#' @importFrom nanny subset
#' 
#' @param .data A tibble
#' @param formula_df A formula dataframe
#' 
make_cens_data = function(.data, formula_df){
	
	# x = enquo(x)
	# alive = enquo(alive)
	
	# For Cibersort
	scale_sd_0_robust = function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y))
	
	# ELIMINATE I HAVE TO SOLVE THIS ISSUE OF NAs
	#####################
	to_eliminate = 
		.data %>% 
		distinct(sample, C, .value_relative) %>%
		group_by(C) %>% 
		mutate(.value_relative = .value_relative %>% boot::logit() %>% scale(scale = F)) %>%
		spread( C, .value_relative) %>%
		as_matrix(rownames = sample) %>%
		{ rownames(.)[apply(., 2, function(x) which(is.na(x))) %>% unlist() %>% as.numeric() %>% unique()] }
	
	.data = 
		.data %>% 
		
		# Eliminate samples with NA
		filter(sample %in% to_eliminate %>% `!`) %>%
		
		filter(!!as.symbol(formula_df$censored_value_column) %>% is.na %>% `!`)
	
	##########################
	
	.data = .data %>% arrange(sample)
	
	

	
	X  = 
		.data %>%
		rename(proportion = .value_relative) %>%
		arrange(sample) %>%
		nest(prop_df = -C) %>%
		
		# Scale
		mutate(prop_df = map(prop_df, ~ mutate(.x, proportion %>% boot::logit() %>% scale(scale = F)))) %>%
		mutate(design = map(prop_df, ~ model.matrix(formula_df$formula_censored_formatted, data=.x) )) %>%
		pull(design) %>%
		abind::abind(along=3)  %>% 
		aperm(perm = c(3, 1, 2)) 
	
	
	S = .data %>% distinct(sample) %>% nrow
	C = .data %>% distinct(C) %>% nrow
	A = ncol(X[1,,])
	
	sample_subset =  .data %>% nanny::subset(sample) %>% arrange(sample)
	
	time = sample_subset %>% mutate(time = !!as.symbol(formula_df$censored_value_column) %>% log1p %>% scale %>% as.numeric) %>% pull(time) %>% as.array()

	#	print(.y)
	which_censored = sample_subset %>% pull(!!as.symbol(formula_df$censored_column)) %>% equals(1) %>% which() %>% as.array()
	which_non_censored = sample_subset %>% pull(!!as.symbol(formula_df$censored_column)) %>% equals(0) %>% which() %>% as.array()
	n_cens = length(which_censored)
	n_non_cens = length(which_non_censored)
	
	list(
		S = S,
		C = C,
		A = A,
		time = time,
		X = X,
		which_censored = which_censored,
		which_non_censored = which_non_censored,
		n_cens = n_cens,
		n_non_cens = n_non_cens,
		prop_C_names = .data %>% distinct(C) %>% arrange(C) %>% pull(C)
	)
	
}

censored_regression = function(.proportions, sampling = F, formula_df, filter_how_many = Inf){
	  
	# x = as.symbol(x)
	# alive = as.symbol(alive)
	
	.proportions %>%
		
		select(sample, formula_df$components_formatted, level, !!formula_df$censored_column, node, C, .draws) %>%
		filter(node %>% is.na %>% `!`) %>%
		unnest(.draws) %>%
		nest(data = -c(node , level,  .chain, .iteration, .draw)) %>%
		
		# Sample half if sampling FALSE
		when(sampling == F ~ (.) %>% group_by(node) %>% sample_frac(0.5) %>% ungroup(), ~ (.)) %>%
		
		when(filter_how_many < nrow(.) ~ (.) %>% sample_n(filter_how_many), ~ (.)) %>%
		# Create input for the model
		mutate(input = imap(data, ~ { make_cens_data(.x, formula_df)})) %>%
		
		# Run model
		mutate(cens_regression = imap(input, ~{ print(.y); run_censored_model(.x, sampling)})) %>%
		select(-data , - input) %>%
		unnest(cens_regression) %>%
		rename(.value = value) %>%
		select(level, node, C, A, .chain, .iteration, .draw, .value, one_of(".draw2", ".lower", ".upper", "prob_non_0")) %>%
		distinct() 
 
}

run_censored_model_joint = function(.data, sampling = F){
	
	sampling_iter = 5
	
	rstan::sampling(
		stanmodels$censored_regression,
		data = .data, chain = 1, iter=150+sampling_iter, warmup=150, save_warmup=F, init="0", refresh = 2000)   %>%
		tidybayes::gather_draws(alpha[A, C]) %>%
		ungroup() %>%
		rename(value = .value, .draw2 = .draw) %>%
		select(-.variable, -.chain, -.iteration)
	
	# %>%
	# 	
	# 	# Change C 
	# 	nest(data = -C) %>%
	# 	mutate(C = !!.data$prop_C_names %>% as.integer)  %>%
	# 	unnest(data) 
	
	
}

#' @keywords internal
#' @importFrom nanny subset
#' 
#' @param .data A tibble
#' @param formula_df A formula dataframe
#' @param relative A boolean
#' 
make_cens_data_joint = function(.data, formula_df, relative = TRUE){
	
	# x = enquo(x)
	# alive = enquo(alive)
	
	# For Cibersort
	scale_sd_0_robust = function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y))
	
	my_value_column = ifelse(relative, as.symbol(".value_relative"), as.symbol(".value_absolute"))
	
	# ELIMINATE I HAVE TO SOLVE THIS ISSUE OF NAs
	#####################
	to_eliminate = 
		.data %>% 
		distinct(sample, new_C, !!my_value_column) %>%
		group_by(new_C) %>% 
		mutate(!!my_value_column := !!my_value_column %>% boot::logit() %>% scale(scale = F)) %>%
		spread( new_C, !!my_value_column) %>%
		as_matrix(rownames = sample) %>%
		{ rownames(.)[apply(., 2, function(x) which(is.na(x))) %>% unlist() %>% as.numeric() %>% unique()] }
	
	.data = 
		.data %>% 
		
		# Eliminate samples with NA
		filter(sample %in% to_eliminate %>% `!`) %>%
		
		filter(!!as.symbol(formula_df$censored_value_column) %>% is.na %>% `!`)
	
	##########################
	
	.data = .data %>% arrange(sample)
	
	
	
	
	X  = 
		.data %>%
		rename(proportion = !!my_value_column) %>%
		arrange(sample) %>%
		nest(prop_df = -new_C) %>%
		
		# Scale
		mutate(prop_df = map(prop_df, ~ mutate(.x, proportion %>% boot::logit() %>% scale(scale = F)))) %>%
		mutate(design = map(prop_df, ~ model.matrix(formula_df$formula_censored_formatted, data=.x) )) %>%
		pull(design) %>%
		abind::abind(along=3)  %>% 
		aperm(perm = c(3, 1, 2)) 
	
	
	S = .data %>% distinct(sample) %>% nrow
	C = .data %>% distinct(new_C) %>% nrow
	A = ncol(X[1,,])
	
	sample_subset =  .data %>% nanny::subset(sample) %>% arrange(sample)
	
	time = sample_subset %>% mutate(time = !!as.symbol(formula_df$censored_value_column) %>% log1p %>% scale %>% as.numeric) %>% pull(time) %>% as.array()
	
	#	print(.y)
	which_censored = sample_subset %>% pull(!!as.symbol(formula_df$censored_column)) %>% equals(1) %>% which() %>% as.array()
	which_non_censored = sample_subset %>% pull(!!as.symbol(formula_df$censored_column)) %>% equals(0) %>% which() %>% as.array()
	n_cens = length(which_censored)
	n_non_cens = length(which_non_censored)
	
	list(
		S = S,
		C = C,
		A = A,
		time = time,
		X = X,
		which_censored = which_censored,
		which_non_censored = which_non_censored,
		n_cens = n_cens,
		n_non_cens = n_non_cens,
		prop_C_names = .data %>% distinct(C) %>% arrange(C) %>% pull(C)
	)
	
}

#' @keywords internal
#' @import parallel
#' 
#' @param .proportions A tibble
#' @param sampling A boolean
#' @param formula_df A formula dataframe
#' @param filter_how_many An integer
#' @param partitions An integer
#' @param relative A boolean
#' 
#' 
censored_regression_joint = 
	function(.proportions, sampling = F, formula_df, filter_how_many = Inf, partitions = 4, relative = TRUE){
	
	# x = as.symbol(x)
	# alive = as.symbol(alive)
	

	partitions = 
		.proportions %>%
		
		# If absolute keep only last level
		#	when(!relative ~ filter(., level == max(level)), ~ (.)) %>%
			
		select(sample, formula_df$components_formatted, level, !!formula_df$censored_column, node, C, .draws) %>%
		filter(node %>% is.na %>% `!`) %>%
		unnest(.draws) %>%
		nest(data = -c(node , level,  .chain, .iteration, .draw)) %>%
		
		# Sample half if sampling FALSE
		when(sampling == F ~ (.) %>% group_by(node) %>% sample_frac(0.5) %>% ungroup(), ~ (.)) %>%
		when(filter_how_many < nrow(.) ~ (.) %>% sample_n(filter_how_many), ~ (.)) %>%
		
		# Add indexes
		rownames_to_column(var = "idx") %>% 
		
		# slice(1:2) %>%
		
		unnest(data) %>%
		unite("idx_C", c(idx, C), remove = F) %>%
		mutate(idx_C = factor(idx_C)) %>%
		#mutate(new_C = as.integer(idx_C)) %>%
		rename(.value_absolute = .value) %>%
		
		# Parallelise
		left_join( (.) %>% distinct(idx_C) %>% mutate(partition = sample(1:partitions, size = n(), replace = T)), by="idx_C") %>%
		nest(data = -partition) 
	
	core_fx = function(.data){
		.data %>%
			nest(data_part = -idx_C) %>%
			mutate(new_C =  1:n()) %>%
			unnest(data_part) %>%
			
			left_join(
				(.) %>%
					# Create input for the model
					make_cens_data_joint(formula_df, relative = relative) %>%
					
					# Run model
					run_censored_model_joint(sampling) %>%
					rename(.value = value) ,
				by=c("new_C" = "C")
			) %>%
			
			select(level, node, C, A, .chain, .iteration, .draw, .value, one_of(".draw2", ".lower", ".upper", "prob_non_0")) %>%
			distinct() 
	}
	
	# Make it parallel
	partitions %>%
		when(
			Sys.info()[['sysname']] != "Windows" ~ (.) %>% pull(data) %>% mclapply(core_fx, mc.cores = 4),
			~ {
				cl <- makeCluster(4)
				(.) %>% pull(data) %>% parlapply(cl, ., core_fx, mc.cores = 4)
			}
		) %>%
		
		# Reformat as tibble
		as_tibble_col(column_name = "data") %>%
		rowid_to_column(var = "partition") %>%
		unnest(data)

}

prepare_TCGA_input = function(file_name, my_dir){
	 
		readRDS(sprintf("%s/TCGA_harmonised/%s", my_dir, file_name)) %>%
		as_tibble(rownames = "ens") %>%
		ensembl_to_symbol(ens) %>%

		gather(sample, count, -ens, -transcript) %>%
		mutate(count = as.integer(count)) %>%
		tidyr::extract(sample, into = "sample", regex = "([a-zA-Z0-9]+-[a-zA-Z0-9]+-[a-zA-Z0-9]+-[a-zA-Z0-9]+)") %>%
		
		# Select primary tumour
		inner_join(
			dir(sprintf("%s/TCGA_harmonised_clinical", my_dir), full.names = T) %>% 
				map_dfr(~ .x %>% readRDS %>% distinct(sample, definition, gender))  %>% 
				#filter(definition == "Primary solid Tumour") %>%
				tidyr::extract(sample, into = "sample", regex = "([a-zA-Z0-9]+-[a-zA-Z0-9]+-[a-zA-Z0-9]+-[a-zA-Z0-9]+)") %>%
				tidyr::extract(sample, into = "patient", regex = "([a-zA-Z0-9]+-[a-zA-Z0-9]+-[a-zA-Z0-9]+)", remove=FALSE) 	,
			by = "sample"
		) %>%
		left_join(
			read_csv("dev/survival_TCGA_curated.csv") %>% 
				#select(bcr_patient_barcode, type, PFI.2, PFI.time.2) %>%
				mutate_each(function(x) ifelse(x == "#N/A", NA, x)) %>%
				type_convert() %>%
				
				# Subset otherwise dataset too big
				select(bcr_patient_barcode, type, PFI.time.2, PFI.2, DSS_cr, DSS.time.cr), 
			by = c("patient" = "bcr_patient_barcode")
		) %>%
		select(-sample) %>%
		mutate_if(is.character, as.factor)
}

draws_to_prob_non_zero = function(.data){
	
	#.value = enquo(.value_col)
	
	.data %>%
		mutate(higher = .value > 0, lower = .value < 0) %>% 
		count(higher) %>%
		spread(higher, n) %>%
		
		# Create column if does not exist
		purrr::when(
			("TRUE" %in% colnames(.) %>% `!`) ~ mutate(., `TRUE` = 0),
			("FALSE" %in% colnames(.) %>% `!`) ~ mutate(., `FALSE` = 0),
			~ (.)
		) %>%
		
		# Smaller probability
		mutate(prob = min(`FALSE`, `TRUE`)/sum(`FALSE`, `TRUE`)) %>% 
		
		# Multiply by 2 and invert 
		mutate(prob = 1 - (prob * 2)) %>%
		
		mutate(prob = ifelse(`FALSE`>`TRUE`, -prob, prob)) %>%
		pull(prob)
}

get_generated_quantities_standalone = function(fit, level, internals){
	
	
	S = internals$Q
	Q = internals$Q
	lv = level
	A = dim(data.frame(internals$X))[2]
	
	mod = switch(
		lv,
		stanmodels$generated_quantities_lv1,
		stanmodels$generated_quantities_lv2,
		stanmodels$generated_quantities_lv3,
		stanmodels$generated_quantities_lv4
	)
	
	fit2 = rstan::gqs(
		mod,
		draws =  as.matrix(fit),
		data = internals$tree_properties
	) 

	
	left_join(
		fit2 %>%
			draws_to_tibble("prop_", "Q", "C") %>%
			mutate(Q = as.integer(Q)) %>%
			mutate(.variable = gsub("_rng", "", .variable)) %>%
			separate(.variable, c("par", "node"), remove = F)  %>%
			select(-par) %>%
			nest(rng_prop = -c(node, C)) %>%
			mutate(C = 1:n()),
		
		fit2 %>%
			draws_to_tibble("mu_", "Q", "C") %>%
			mutate(Q = as.integer(Q)) %>%
			mutate(.variable = gsub("_rng", "", .variable)) %>%
			separate(.variable, c("par", "node"), remove = F)  %>%
			select(-par) %>%
			nest(rng_mu = -c(node, C)) %>%
			mutate(C = 1:n()),
		by=c("C", "node")
	)
	
	
}

#' @importFrom nanny permute_nest
#' 
#' @param .data A tibble
#' 
get_signatures = function(.data){
	.data$proportions %>%
		filter(.variable %>% is.na %>% `!`) %>%
		select(-proportions, -rng) %>%
		unnest(draws) %>%
		
		# Group
		nest(node = -c(level, .variable, A)) %>%
		mutate(node = map(
			node,
			~ .x %>%
				
				# Build combination of cell types
				nanny::permute_nest(
					.names_from = `Cell type category`,
					.values_from = c(.value)
				) %>% 
				
				# Perform calculation
				mutate(prob_df = map(
					data, 
					~.x %>%
						nest(data = -c(`Cell type category`)) %>% 
						mutate(med = map_dbl(data, ~.x$.value %>% median )) %>% 
						mutate(med = rev(med)) %>% 
						mutate(frac_up = map2_dbl(data, med, ~ (.x$.value  > .y) %>% sum %>% magrittr::divide_by(nrow(.x)))) %>% 
						mutate(frac_down = map2_dbl(data, med, ~ (.x$.value  < .y) %>% sum %>% magrittr::divide_by(nrow(.x)))) %>%
						mutate(prob = ifelse(frac_up > frac_down, frac_up, frac_down)) %>%
						
						# stretch to o 1 interval
						mutate(prob = (prob-0.5)/0.5) %>%
						
						# Insert sign
						mutate(prob = ifelse(frac_up < frac_down, -prob, prob)) %>%
						select(`Cell type category`, prob)
				)) %>%
				select(-data) %>%
				unnest(prob_df) %>% 
				filter(`Cell type category_1` == `Cell type category`) %>%
				select(-`Cell type category`)
			
		)) %>%
		unnest( node)
}

get_CI = function(.data, credible_interval = 0.90, cluster_CI = 0.55) {
	
	# Choose the test
	if("draws_cens" %in% colnames(.data$proportions)){
		.v = as.symbol(".value_2")
		.l = as.symbol(".lower_2")
		.u = as.symbol(".upper_2")
	} else {
		.v = as.symbol(".value_alpha2")
		.l = as.symbol(".lower_alpha2")
		.u = as.symbol(".upper_alpha2")
	}
	
	.d = 
		.data$proportions %>%
		filter(.variable %>% is.na %>% `!`) %>%
		cluster_posterior_slopes(credible_interval = cluster_CI) %>%
		extract_CI(credible_interval)
	
	dx = list()	
	
	# Level 1
	if(.d %>% filter(level ==1) %>% nrow %>% `>` (0))
		dx = 
		.d %>%
		
		filter(level ==1) %>%
		mutate(fold_change_ancestor = 0) %>%
		identify_baseline_by_clustering( ) %>%
		
		mutate(significant = ((!!.l - 0) * (!!.u - 0)) > 0) %>%
		mutate(fold_change  = ifelse(significant, !!.v, 0))
	
	# Level 2
	if(.d %>% filter(level ==2) %>% nrow %>% `>` (0))
		dx =	dx %>% bind_rows(
			.d %>%
				
				filter(level ==2) %>%
				left_join(ancestor_child, by = "Cell type category") %>%
				left_join( dx %>%	select(ancestor = `Cell type category`, fold_change_ancestor = fold_change),  by = "ancestor" )  %>%
				identify_baseline_by_clustering( ) %>%
				
				mutate(significant = ((!!.l - 0) * (!!.u - 0)) > 0) %>%
				mutate(fold_change  = ifelse(significant, !!.v, 0))
		) 
	
	
	# Level 3
	if(.d %>% filter(level ==3) %>% nrow %>% `>` (0))
		dx =	dx %>% bind_rows(
			.d %>%
				
				filter(level ==3) %>%
				left_join(ancestor_child, by = "Cell type category") %>%
				left_join( dx %>%	select(ancestor = `Cell type category`, fold_change_ancestor = fold_change) ,  by = "ancestor")  %>%
				identify_baseline_by_clustering( ) %>%
				
				mutate(significant = ((!!.l - 0) * (!!.u - 0)) > 0) %>%
				mutate(fold_change  = ifelse(significant, !!.v, 0))
		)
	
	# Level 4
	if(.d %>% filter(level ==4) %>% nrow %>% `>` (0))
		dx =	dx %>% bind_rows(
			.d %>%
				
				filter(level ==4) %>%
				left_join(ancestor_child, by = "Cell type category") %>%
				left_join( dx %>%	select(ancestor = `Cell type category`, fold_change_ancestor = fold_change) ,  by = "ancestor")  %>%
				identify_baseline_by_clustering( ) %>%
				
				mutate(significant = ((!!.l - 0) * (!!.u - 0)) > 0) %>%
				mutate(fold_change  = ifelse(significant, !!.v, 0))
		)
	
	dx
	
	
}

