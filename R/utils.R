


# library(magrittr)
# library(tidyverse)
# library(foreach)
# library(rstan)
# library(tidyTranscriptomics)

# gtools::rdirichlet(150, c(200, 200, 1, 1)) %>% sirt::dirichlet.mle()

#' This is a generalisation of ifelse that acceots an object and return an objects
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

#' format_for_MPI
#'
#' @description Format reference data frame for MPI
format_for_MPI = function(df, shards) {
	df %>%

		left_join((.) %>%
								distinct(G) %>%
								arrange(G) %>%
								mutate(idx_MPI = head(
									rep(1:shards, (.) %>% nrow %>% `/` (shards) %>% ceiling), n = (.) %>% nrow
								))) %>%
		arrange(idx_MPI, G) %>%

		# Decide start - end location
		group_by(idx_MPI) %>%
		do((.) %>%
			 	left_join(
			 		(.) %>%
			 			distinct(sample, G) %>%
			 			arrange(G) %>%
			 			count(G) %>%
			 			mutate(end = cumsum(n)) %>%
			 			mutate(start = c(
			 				1, .$end %>% rev() %>% `[` (-1) %>% rev %>% `+` (1)
			 			))
			 	)) %>%
		ungroup() %>%

		# Add ct_symbol MPI rows indexes - otherwise spread below gives error
		left_join(
			(.) %>%
				group_by(idx_MPI) %>%
				distinct(G) %>%
				arrange(G) %>%
				mutate(`symbol MPI row` = 1:n()) %>%
				ungroup
		) %>%

		# Add counts MPI rows indexes
		group_by(idx_MPI) %>%
		arrange(G) %>%
		mutate(`read count MPI row` = 1:n()) %>%
		# do( (.) %>% arrange(G) %>% rowid_to_column("read count MPI row") ) %>%
		ungroup

}

#' xx
#' @import magrittr
format_for_MPI_from_linear = function(df) {
	shards = df %>% arrange(shards %>% desc) %>% slice(1) %>% pull(shards)
	if (shards %>% length %>% equals(0))
		shards = 1


	df %>%

		left_join((.) %>%
								distinct(GM) %>%
								arrange(GM) %>%
								mutate(idx_MPI = head(
									rep(1:shards, (.) %>% nrow %>% `/` (shards) %>% ceiling), n = (.) %>% nrow
								))) %>%
		arrange(idx_MPI, GM, G) %>%

		# Add counts MPI rows indexes
		group_by(idx_MPI) %>%
		arrange(GM, G) %>%
		do((.) %>% rowid_to_column("read count MPI row")) %>%
		ungroup

}

format_for_MPI_from_linear_GM = function(df) {
	shards = df %>% arrange(shards %>% desc) %>% slice(1) %>% pull(shards)
	if (shards %>% length %>% equals(0))
		shards = 1

	df %>%

		left_join((.) %>%
								distinct(GM) %>%
								arrange(GM) %>%
								mutate(idx_MPI = head(
									rep(1:shards, (.) %>% nrow %>% `/` (shards) %>% ceiling), n = (.) %>% nrow
								))) %>%
		arrange(idx_MPI, GM) %>%

		# Add counts MPI rows indexes
		group_by(idx_MPI) %>%
		arrange(GM) %>%
		do((.) %>% rowid_to_column("read count MPI row")) %>%
		ungroup

}

format_for_MPI_from_linear_dec = function(df, lv) {
	shards = df %>% arrange(shards %>% desc) %>% slice(1) %>% pull(shards)
	if (shards %>% length %>% equals(0))
		shards = 1

	df %>%

		left_join((.) %>%
								distinct(GM) %>%
								arrange(GM) %>%
								mutate(idx_MPI = head(
									rep(1:shards, (.) %>% nrow %>% `/` (shards) %>% ceiling), n = (.) %>% nrow
								))) %>%
		arrange(idx_MPI, GM, !!as.symbol(sprintf("C%s", lv))) %>%

		# Add counts MPI rows indexes
		group_by(idx_MPI) %>%
		arrange(GM, !!as.symbol(sprintf("C%s", lv))) %>%
		do((.) %>% rowid_to_column("read count MPI row")) %>%
		ungroup

}

#' add_partition
#'
#' @description Add partition column dto data frame
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

#' get_MPI_deconv
#'
#' @description Get data format for MPI deconvolution part
get_MPI_deconv = function(y_source, shards, my_level, tree) {
	# This function is needed in case
	# I have too many shards and not enough data
	add_empty_shards = function(df) {
		tibble(partition = 1:shards, n = 0 %>% as.integer) %>%
			anti_join(df, by = "partition") %>%
			bind_rows(df) %>%
			arrange(partition)
	}

	y_MPI_source =
		y_source %>%
		distinct(level, Q, S, symbol, G, GM, `Cell type category`, `read count`) %>%
		filter(level == my_level)  %>%

		# Add universal cell type rank
		left_join(tree %>%
								data.tree::ToDataFrameTree("name", sprintf("C%s", my_level)) %>%
								select(-1) %>%
								setNames(c("Cell type category" , "ct_rank"))) %>%

		# Arrange very important for consistency
		arrange(Q, symbol, ct_rank) %>%
		add_partition("symbol", shards) %>%
		mutate(partition = partition %>% as.integer) %>%
		group_by(partition) %>%
		left_join((.) %>% distinct(symbol) %>% mutate(MPI_row = 1:n())) %>%
		ungroup()

	list(
		y_MPI_source = y_MPI_source,

		y_MPI_symbol_per_shard =
			y_MPI_source %>%
			distinct(symbol, partition) %>%
			count(partition) %>%
			add_empty_shards %>%
			spread(partition, n) %>%
			as_vector %>% array,

		y_MPI_idx_symbol =
			y_MPI_source %>%
			distinct(MPI_row, GM, partition, Q) %>%
			spread(partition, GM) %>%
			select(-Q) %>% distinct %>%
			replace(is.na(.), 0 %>% as.integer) %>%
			as_matrix(rownames = "MPI_row") %>%
			t,

		y_MPI_G_per_shard =
			y_MPI_source %>%
			distinct(symbol, `Cell type category`, partition) %>%
			count(partition) %>%
			add_empty_shards %>%
			spread(partition, n) %>%
			as_vector %>% array,

		y_MPI_idx =
			y_MPI_source %>%
			distinct(partition, symbol, G, `Cell type category`, MPI_row) %>%
			select(-symbol) %>%
			spread(partition, G) %>%
			arrange(MPI_row, `Cell type category`) %>%
			select(-MPI_row,-`Cell type category`) %>%
			replace(is.na(.), 0 %>% as.integer) %>%
			as_matrix %>%
			t,

		y_idx =
			y_MPI_source %>% distinct(symbol, G, `Cell type category`) %>% arrange(symbol, `Cell type category`) %>% pull(G),

		y_MPI_N_per_shard =
			y_MPI_source %>%
			distinct(MPI_row, `read count`, partition, Q) %>%
			count(partition) %>%
			add_empty_shards %>%
			spread(partition, n) %>%
			as_vector %>% array,

		y_MPI_count =
			y_MPI_source %>%
			distinct(MPI_row, `read count`, partition, Q) %>%
			spread(partition, `read count`) %>%
			arrange(MPI_row, Q) %>%
			select(-MPI_row,-Q) %>%
			replace(is.na(.), 0 %>% as.integer) %>%
			as_matrix %>%
			t
	)
}

#' Plot differences between inferred and observed transcription abundances
#'
#' @description  Plot differences between inferred and observed transcription abundances
plot_differences_in_lambda = function() {
	# Plot differences in lambda_log
	(
		fit %>%
			tidybayes::gather_draws(lambda_log[G]) %>%
			tidybayes::median_qi() %>%
			left_join(
				counts_baseline %>%
					distinct(`symbol`, G, `Cell type category`)
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
	)  %>% plotly::ggplotly()

	#
	(
		fit %>%
			extract(
				pars = c(
					"lambda_mu",
					"lambda_sigma",
					"exposure_rate",
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
	) %>% plotly::ggplotly()

}

#' Print which reference genes are in the mix
#'
#' @description Print which reference genes are in the mix
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

#' plot_counts_inferred_sum
#'
#' @description Get data format for MPI deconvolution part
plot_counts_inferred_sum = function(fit_obj, samples = NULL) {
	fit_obj %$% fit %>%
		rstan::summary(par = c("nb_sum")) %$% summary %>%
		as_tibble(rownames = "par") %>% select(par, `2.5%`, `50%`, `97.5%`) %>%
		separate(par, c(".variable", "Q", "GM"), sep = "\\[|,|\\]") %>%
		mutate(Q = Q %>% as.integer, GM = GM %>% as.integer) %>%
		left_join(fit_obj %$% data_source %>% distinct(Q, GM, symbol, `read count`, sample)) %>%

		# Select samples
		{
			if (samples %>% is.null %>% `!`)
				(.) %>% filter(sample %in% !samples)
			else
				(.)
		} %>%

		# Check if inside
		rowwise %>%
		mutate(inside = between(`read count`, `2.5%`, `97.5%`)) %>%
		ungroup %>%
		ggplot(aes(
			x = `read count` + 1,
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

#' choose_chains_majority_roule
#'
#' @description Get which chain cluster is more opulated in case I have divergence
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

#' filter_reference
#'
#' @description Filter the reference
#'
#' @export
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
		filter(symbol %in% (mix %>% colnames)) %>%
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
							ungroup
					) %>%

					# Print number of markewrs per comparison
					{
						(.) %>% distinct(symbol, ct1, ct2) %>% count(ct1, ct2) %>% print(n = 99)
						(.)
					},

				# Get house keeping genes
				(.) %>% filter(`house keeping`)
			)
		}
}

#' get_idx_level
#'
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

#' parse_summary
#'
#' @description Parse the stan fit object
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

#' parse_summary_check_divergence
#'
#' @description Parse the stan fit object and check for divergencies
parse_summary_check_divergence = function(draws) {
	draws %>%

		group_by(level, .variable,  Q, sample,  C, `Cell type category`) %>%

		# If not converged choose the majority chains
		mutate(converged = diptest::dip.test(`.value_relative`) %$%	`p.value` > 0.05) %>%

		# If some proportions have not converged chose the most populated one
		do(
			(.) %>%
				ifelse_pipe(
					(.) %>% distinct(converged) %>% pull(1) %>% `!`,
					~ .x %>% choose_chains_majority_roule
				)
		) %>%

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

#' create_tree_object
#'
#' @description create tree object that is in data directory
#'
#' @export
create_tree_object = function(my_ref = ARMET::ARMET_ref) {
	#yaml:: yaml.load_file("~/PhD/deconvolution/ARMET/data/tree.yaml") %>%
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
}

#' as_matrix
#'
#' @description Convert tibble to matrix
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

#' ToDataFrameTypeColFull
#'
#' @description Extension of data.tree package. It converts the tree into data frame
#'
#' @export
ToDataFrameTypeColFull = function(tree, fill = T, ...) {
	tree %>%
		data.tree::Clone() %>%
		{
			t = (.)
			foreach(l = 1:(t %$% Get("level") %>% max), .combine = bind_rows) %do% {
				data.tree::Clone(t) %>%
					{
						data.tree::Prune(., function(x)
							x$level <= l + 1)
						.
					} %>%
					data.tree::ToDataFrameTypeCol(...) %>%
					as_tibble
			}
		} %>%
		distinct() %>%
		ifelse_pipe(
			fill,
			~ .x %>%
				{
					if ("level_3" %in% ((.) %>% colnames))
						(.) %>% mutate(level_3 = ifelse(level_3 %>% is.na, level_2, level_3))
					else
						(.)
				} %>%
				{
					if ("level_4" %in% ((.) %>% colnames))
						(.) %>% mutate(level_4 = ifelse(level_4 %>% is.na, level_3, level_4))
					else
						(.)
				} %>%
				{
					if ("level_5" %in% ((.) %>% colnames))
						(.) %>% mutate(level_5 = ifelse(level_5 %>% is.na, level_4, level_5))
					else
						(.)
				} %>%
				{
					if ("level_6" %in% ((.) %>% colnames))
						(.) %>% mutate(level_6 = ifelse(level_6 %>% is.na, level_5, level_6))
					else
						(.)
				}
			) %>%
		select(..., everything())
}

#' vb_iterative
#'
#' @description Runs iteratively variational bayes until it suceeds
#'
#' @importFrom rstan vb
#'
#' @param model A Stan model
#' @param output_samples An integer of how many samples from posteriors
#' @param iter An integer of how many max iterations
#' @param tol_rel_obj A real
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

parse_baseline = function(.data, shards_in_levels, lv) {
	.data %>%
		filter(level == lv) %>%
		distinct(
			sample,
			symbol,
			`Cell type category`,
			level,
			`read count`,
			counts_idx,
			G,
			GM,
			S,
			`house keeping`
		) %>%
		left_join(tibble(level = lv, shards = shards_in_levels)) %>%
		format_for_MPI_from_linear()
}

get_MPI_df = function(counts_baseline_to_linear,
											y_source,
											counts_baseline,
											shards_in_levels,
											lv) {
	list(
		counts_idx_lv_MPI =
			counts_baseline_to_linear %>%
			parse_baseline(shards_in_levels, lv)  %>%
			distinct(idx_MPI, counts_idx, `read count MPI row`) %>%
			spread(idx_MPI,  counts_idx) %>%
			select(-`read count MPI row`) %>%
			replace(is.na(.), -999 %>% as.integer) %>%
			as_matrix() %>% t %>% 		as.data.frame,

		size_counts_idx_lv_MPI =
			counts_baseline_to_linear %>%
			parse_baseline(shards_in_levels, lv)   %>%
			distinct(idx_MPI, counts_idx, `read count MPI row`) %>%
			count(idx_MPI) %>%
			pull(n) %>%
			ifelse_pipe(length((.)) == 0, ~ 0) %>%  		as.array,

		# Count indexes
		counts_G_lv_MPI =
			counts_baseline_to_linear %>%
			parse_baseline(shards_in_levels, lv)   %>%
			distinct(idx_MPI, G, `read count MPI row`)  %>%
			spread(idx_MPI,  G) %>%
			select(-`read count MPI row`) %>%
			replace(is.na(.), -999 %>% as.integer) %>%
			as_matrix() %>% t %>% 		as.data.frame,

		size_counts_G_lv_MPI =
			counts_baseline_to_linear %>%
			parse_baseline(shards_in_levels, lv)   %>%
			distinct(idx_MPI, G, `read count MPI row`)  %>%
			count(idx_MPI) %>%
			pull(n) %>%
			ifelse_pipe(length((.)) == 0, ~ 0) %>%  		as.array,

		counts_G_lv_MPI_non_redundant =
			counts_baseline_to_linear %>%
			parse_baseline(shards_in_levels, lv)   %>%
			distinct(idx_MPI, G)  %>%
			group_by(idx_MPI) %>% do((.) %>% rowid_to_column("read count MPI row")) %>% ungroup() %>%
			spread(idx_MPI,  G) %>%
			select(-`read count MPI row`) %>%
			replace(is.na(.), -999 %>% as.integer) %>%
			as_matrix() %>% t %>% 		as.data.frame,

		size_counts_G_lv_MPI_non_redundant =
			counts_baseline_to_linear %>%
			parse_baseline(shards_in_levels, lv)   %>%
			distinct(idx_MPI, G)  %>%
			count(idx_MPI) %>%
			pull(n) %>%
			ifelse_pipe(length((.)) == 0, ~ 0) %>%  		as.array,

		counts_G_lv_MPI_non_redundant_reps =
			counts_baseline_to_linear %>%
			parse_baseline(shards_in_levels, lv)   %>%
			distinct(idx_MPI, G, `read count MPI row`)  %>%
			left_join((.) %>% count(idx_MPI, G)) %>%
			distinct(idx_MPI, G, n) %>%
			group_by(idx_MPI) %>% do((.) %>% rowid_to_column("read count MPI row")) %>% ungroup() %>%
			distinct(idx_MPI, n, `read count MPI row`) %>%
			spread(idx_MPI,  n) %>%
			select(-`read count MPI row`) %>%
			replace(is.na(.), -999 %>% as.integer) %>%
			as_matrix() %>% t %>% 		as.data.frame,

		# Count indexes
		counts_S_lv_MPI =
			counts_baseline_to_linear %>%
			parse_baseline(shards_in_levels, lv) %>%
			distinct(idx_MPI, S, `read count MPI row`)  %>%
			spread(idx_MPI,  S) %>%
			select(-`read count MPI row`) %>%
			replace(is.na(.), -999 %>% as.integer) %>%
			as_matrix() %>% t %>% 		as.data.frame,

		size_counts_S_lv_MPI =
			counts_baseline_to_linear %>%
			parse_baseline(shards_in_levels, lv) %>%
			distinct(idx_MPI, S, `read count MPI row`)   %>%
			count(idx_MPI) %>%
			pull(n) %>%
			ifelse_pipe(length((.)) == 0, ~ 0) %>%  		as.array,

		# mix Counts
		y_linear_MPI =
			y_source %>%
			filter(level == lv) %>%
			distinct(GM, Q, S, `read count`, level) %>%
			left_join(tibble(level = lv, shards = shards_in_levels)) %>%
			format_for_MPI_from_linear_GM() %>%
			distinct(idx_MPI, `read count`, `read count MPI row`)  %>%
			spread(idx_MPI,  `read count`) %>%
			select(-`read count MPI row`) %>%
			replace(is.na(.), -999 %>% as.integer) %>%
			as_matrix() %>% t %>% 		as.data.frame,

		size_y_linear_MPI =
			y_source %>%
			filter(level == lv) %>%
			distinct(GM, Q, S, `read count`, level) %>%
			left_join(tibble(level = lv, shards = shards_in_levels)) %>%
			format_for_MPI_from_linear_GM() %>%
			distinct(idx_MPI, `read count`, `read count MPI row`)  %>%
			count(idx_MPI) %>%
			pull(n) %>%
			ifelse_pipe(length((.)) == 0, ~ 0) %>%  		as.array,

		y_linear_S_MPI =
			y_source %>%
			filter(level == lv) %>%
			distinct(GM, Q, S, `read count`, level) %>%
			left_join(tibble(level = lv, shards = shards_in_levels)) %>%
			format_for_MPI_from_linear_GM() %>%
			distinct(idx_MPI, S, `read count MPI row`)  %>%
			spread(idx_MPI,  S) %>%
			select(-`read count MPI row`) %>%
			replace(is.na(.), -999 %>% as.integer) %>%
			as_matrix() %>% t %>% 		as.data.frame,

		size_y_linear_S_MPI =
			y_source %>%
			filter(level == lv) %>%
			distinct(GM, Q, S, `read count`, level) %>%
			left_join(tibble(level = lv, shards = shards_in_levels)) %>%
			format_for_MPI_from_linear_GM() %>%
			distinct(idx_MPI, S, `read count MPI row`)  %>%
			count(idx_MPI) %>%
			pull(n) %>%
			ifelse_pipe(length((.)) == 0, ~ 0) %>%  		as.array,

		G_linear_MPI =
			counts_baseline %>% filter(level == lv) %>%

			# I have fixed this for right order
			select(level, G, GM, sprintf("C%s", lv)) %>%
			distinct() %>%
			arrange(GM, !!as.symbol(sprintf("C%s", lv))) %>%

			#distinct(G, GM, C, level) %>%
			left_join(tibble(level = lv, shards = shards_in_levels)) %>%
			format_for_MPI_from_linear_dec(lv) %>%
			distinct(idx_MPI, G, `read count MPI row`) %>%
			spread(idx_MPI,  G) %>%
			select(-`read count MPI row`) %>%
			replace(is.na(.), -999 %>% as.integer) %>%
			as_matrix() %>% t %>% 		as.data.frame,

		size_G_linear_MPI =
			counts_baseline %>% filter(level == lv) %>%
			# I have fixed this for right order
			select(level, G, GM, sprintf("C%s", lv)) %>%
			distinct() %>%
			arrange(GM, !!as.symbol(sprintf("C%s", lv))) %>%

			#distinct(G, GM, C, level) %>%
			left_join(tibble(level = lv, shards = shards_in_levels)) %>%
			format_for_MPI_from_linear_dec(lv) %>%
			distinct(idx_MPI, G, `read count MPI row`)   %>%
			count(idx_MPI) %>%
			pull(n) %>%
			ifelse_pipe(length((.)) == 0, ~ 0) %>%  		as.array
	)
}

ref_mix_format = function(ref, mix) {
	bind_rows(
		# Get reference based on mix genes
		ref %>% mutate(`query` = FALSE),
		mix %>%
			gather(`symbol`, `read count`,-sample) %>%
			inner_join(ref %>% distinct(symbol)) %>%
			left_join(ref %>% distinct(symbol, `house keeping`)) %>%
			mutate(`Cell type category` = "query") %>%
			mutate(`query` = TRUE)
	)	%>%

		# Add marker symbol indeces
		left_join((.) %>%
								filter(!`house keeping`) %>%
								distinct(`symbol`) %>%
								mutate(M = 1:n())) %>%

		# Add sample indeces
		arrange(!`query`) %>% # query first
		mutate(S = factor(sample, levels = .$sample %>% unique) %>% as.integer) %>%

		# Add query samples indeces
		left_join((.) %>%
								filter(`query`) %>%
								distinct(`sample`) %>%
								mutate(Q = 1:n())) %>%

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
				ungroup()
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
				mutate(G = 1:n())
		) %>%
		left_join(
			(.) %>%
				filter(!`house keeping` & !`query`) %>%
				distinct(level, symbol) %>%
				arrange(level, symbol) %>%
				mutate(GM = 1:n()) %>%
				select(-level)
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

#' @export
draws_to_alphas = function(.data, pars) {
	.data %>%
		tidybayes::gather_draws(`prop_[1,a-z]`[Q, C], regex = T) %>%
		ungroup() %>%
		{
			print((.))
			(.)
		} %>%
		filter(.variable %in% pars) %>%
		{
			print((.))
			(.)
		} %>%
		nest(data = -c(.variable, Q)) %>%
		mutate(
			alphas = map(
				data,
				~
					.x %>%
					spread(C, .value) %>%
					select(-c(1:3)) %>%
					as_matrix() %>%
					sirt::dirichlet.mle() %$%
					alpha %>%
					as_tibble() %>%
					rename(alpha = value) %>%
					mutate(C = 1:n())
			)
		) %>%
		select(-data) %>%
		unnest(cols = alphas) %>%
		spread(C, alpha) %>%
		select(-Q) %>%
		nest(alphas = -.variable) %>%
		pull(alphas)
}

draws_to_exposure = function(.data) {
	.data %>%
		tidybayes::gather_draws(exposure_rate[Q]) %>%
		summarise(.mean = .value %>% mean, .sd = .value %>% sd) %>%
		ungroup() %>%
		select(.mean , .sd)
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

plot_markers = function(.data, n_markers, mix, ct1, ct2, n = 10, level) {
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

#' This is a generalisation of ifelse that acceots an object and return an objects
#'
#' @import ggplot2
#'
#' @export
level_to_plot_inferred_vs_observed  = function(result, level, S = NULL, cores = 20){

	library(multidplyr)

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
#'
#' @param input.df A tibble
#' @param formula A formula
#' @param sample_column A symbol
#' @export
create_design_matrix = function(input.df, formula, sample_column){

	sample_column = enquo(sample_column)

	model.matrix(
		object = formula,
		data =
			input.df %>%
			select(!!sample_column, one_of(parse_formula(formula))) %>%
			distinct %>% arrange(!!sample_column)

	)

}

#' Formula parser
#'
#' @param fm A formula
#'
#' @return A character vector
#'
#'
parse_formula <- function(fm) {
	if (attr(terms(fm), "response") == 1)
		stop("The formula must be of the kind \"~ covariates\" ")
	else
		as.character(attr(terms(fm), "variables"))[-1]
}

rebuild_last_component_sum_to_zero = function(.){
	(.) %>%
		group_by(.variable) %>%
		do({
			max_c = (.) %>% pull(C) %>% max
			bind_rows(
				(.) %>% filter(C < max_c | A > 1),
				(.) %>%
					filter(C < max_c & A == 1) %>%
					group_by(.chain, .iteration, .draw , A, .variable ) %>%
					summarise(.value = -sum(.value)) %>%
					mutate(C = max_c)
			)
		}) %>%
		ungroup()
}

get_relative_zero = function(fit_parsed){
	# fit %>%
	# 	tidybayes::gather_draws(`alpha_[1abcdefghi]`[A, C], regex = T) %>%
	# 	ungroup() %>%
	fit_parsed %>%
		#filter(A == 2) %>%
		nest(data = -.variable) %>%
		mutate(zero = map(
			data,
			~ .x %>%
				filter(A == 2) %>%
				pull(.value) %>%
				density(na.rm = T) %>%
				{	tibble(x = (.)$x, y = (.)$y)	} %>%
				arrange(y %>% desc) %>% slice(1) %>% select(zero = x)
		)) %>%
		unnest(cols = c(data, zero))
}
