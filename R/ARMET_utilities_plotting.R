#' Plots the densities distributions columnwise of a matrix
#'
#' @param df A matrix
#' @param color A char
#' @param fill A array of numbers
#' @param alpha A number
#' @param do_log A bool
#' @return A ggplot
plot_densities = function(df, color="0",  fill = "0", alpha = 0.30, do_log = T){

	if(do_log) df = df %>% dplyr::mutate(value=value+0.1)
	p =
		ggplot2::ggplot(	df, ggplot2::aes(value, group=sample, color = factor(color))) +
		ggplot2::geom_line(stat="density", alpha=alpha) +
		ggplot2::expand_limits(x=0.1) +
		ggplot2::theme_bw() +
		ggplot2::theme(
			panel.border = ggplot2::element_blank(),
			panel.grid.major = ggplot2::element_blank(),
			panel.grid.minor = ggplot2::element_blank(),
			axis.line = ggplot2::element_line(colour = "black")
		)

	if(do_log) p = p + ggplot2::scale_x_log10()

	return( p )
}

#' Plots the densities from a melted data frame
#'
#' @param df_melted A matrix
#' @return A ggplot
plot_densities_from_melted = function(df_melted){

	p = ggplot2::ggplot(df_melted, ggplot2::aes(value+0.1, group=X2, color=source)) +
		ggplot2::geom_line(stat="density", alpha=0.15) +
		ggplot2::scale_x_log10() +
		ggplot2::expand_limits(x=0.1) +
		ggplot2::theme_bw() +
		ggplot2::theme(
			panel.border = ggplot2::element_blank(),
			panel.grid.major = ggplot2::element_blank(),
			panel.grid.minor = ggplot2::element_blank(),
			axis.line = ggplot2::element_line(colour = "black")
		)
	return( p )
}

#' Parse the annotated tree
#' @rdname ARMET_plotFit
#'
#' Prints a report of the hipothesis testing
#'
#' @param ARMET-tc object
#' @param cell type string
#'
#' @return a grid pbject
#'
#' @examples
#'  ARMET_plotFit(ARMET_tc_result, ct = "TME")
#' @export
ARMET_plotFit = function(obj, ct = "TME", param = "estimate_prop_with_uncertanties", dodge = 0){

	node_info = get_node_from_name(obj$tree, ct)

	is_categorical = all(
		obj$input$my_design %>%
			tibble::as_tibble() %>%
			dplyr::select(!!obj$input$cov_to_test) %>%
			dplyr::pull() %>% unique() %>% sort() == c(0,1)
	)

	my_noise = rnorm( nrow(node_info[[param]]), 0.1, dodge	)

	my_method = ifelse(is_categorical, "lm", "loess")

	plot_props = ggplot2::ggplot(
		node_info[[param]], ggplot2::aes(x=get(obj$input$cov_to_test)+my_noise, y=estimate, colour=ct, group=ct)) +
		ggplot2::geom_point() +
		ggplot2::geom_errorbar(
			ggplot2::aes(ymin=conf.low, ymax=conf.high, colour=ct),
			width=0, alpha=0.4
		) +
		#ggplot2::scale_colour_brewer(palette = "Set1") +
		#ggplot2::scale_fill_brewer(palette = "Set1") +
		ggplot2::geom_smooth(method = my_method, alpha = 0.05) +
		ggplot2::theme_bw() +
		scale_y_continuous(limits = c(1e-10, 1), trans ='log10' ) +
		ggplot2::ggtitle("Scatter plot of cell type proportion") +
		labs(x = obj$input$cov_to_test, color="Cell type")

	plot_coef_ang = ggplot2::ggplot(
		node_info$coef_ang_posterior_adj , ggplot2::aes(value, group=ct, color=ct)) +
		ggplot2::geom_line(stat="density", size=1) +
		#ggplot2::scale_colour_brewer(palette = "Set1") +
		#ggplot2::scale_fill_brewer(palette = "Set1") +
		ggplot2::expand_limits(x=0.1) +
		ggplot2::theme_bw() +
		ggplot2::theme(
			panel.border = ggplot2::element_blank(),
			axis.line = ggplot2::element_line(colour = "black")
		) +
		ggplot2::ggtitle("Posterior distribution of angular coeff.") +
		labs(x = "Trend of change", color="Cell type")


	gridExtra::grid.arrange(plot_props, plot_coef_ang, nrow=2)

}

#' Create tree plot of results
#' @rdname ARMET_plotTree
#'
#' Prints a report of the hipothesis testing
#'
#' @param ARMET-tc object
#'
#' @return a ggplot
#'
#' @export
ARMET_plotTree = function(obj){

	library(data.tree)
	library(ggtree)
	library(foreach)
	library(ggplot2)
	library(tibble)
	library(dplyr)

	# Build phylo tree
	tt = get_tree_hypoth_test(treeYaml, obj$tree)
	tt_phylo = tt %>% as.phylo

	# Formatted names
	ct_names_formatted =
		ToDataFrameTree(tt, "name") %>%
		as_tibble %>% select(-levelName) %>%
		bind_cols(taxa = c(
			"Root",	"Epi", "Endo", "Fibro", "Adipo", "Immune", "Granulo", "Eosin", "Neutro", "Mono deriv", "Mono", "M0 macro", "M1 macro", "M2 macro", "Dendr rest", "Drondr activ", "Mast", "Activ", "Rest", "B", "Naive","Mem", "T", "CD8", "CD4 naive", "H1", "H2", "Follic", "γδ", "Reg", "Mem cent", "Mem activ", "NK", "Activ", "Rest", "Plasma"
		))

	# Give formatted names
	tt_phylo$tip.label = ct_names_formatted %>% filter(name %in% tt_phylo$tip.label) %>% arrange(match(name, tt_phylo$tip.label)) %>% pull(taxa)
	tt_phylo$node.label = ct_names_formatted %>% filter(name %in% tt_phylo$node.label) %>% arrange(match(name, tt_phylo$node.label)) %>% pull(taxa)

	# Create annotation
	internal_branch_length = 40
	external_branch_length = 10
	dd =
		ToDataFrameTree(
			tt,
			"name",
			"isLeaf", "level",
			"estimate_extrinsic",
			"std_error_extrinsic",
			"direction_extrinsic",
			"pvalue_extrinsic",
			"significance_extrinsic",
			"Cell type proportion"
		) %>%
		as_tibble() %>%
		left_join(ct_names_formatted, by= "name") %>%
		select(-levelName) %>% select(taxa, everything()) %>%
		mutate(estimate_extrinsic = ifelse(pvalue_extrinsic<0.05, estimate_extrinsic, NA)) %>%
		mutate(branch_length = ifelse(isLeaf, 0.1, 2)) %>%

		# Correct branch length
		mutate(level = level -1) %>%
		mutate(branch_length = ifelse(!isLeaf, internal_branch_length,	external_branch_length) ) %>%
		mutate(branch_length =  ifelse(	isLeaf, branch_length + ((max(level) - level) * internal_branch_length),	branch_length	))

	# Change branch length
	tt_phylo$edge.length = (dd %>% arrange(match(taxa, c(tt_phylo$tip.label, tt_phylo$node.label))) %>% pull(branch_length)) [tt_phylo$edge[,2]]

	# Plot square tree
	# tt_phylo %>%
	# 	ggtree(aes(color=estimate_extrinsic, size=`Cell type proportion`)) %<+% dd +
	# 	xlim(-10, 150) +
	# 	geom_tiplab(size=5, hjust = 0, color="black" , offset = 2) +
	# 	geom_nodelab(size=4, hjust = 1, nudge_x = -3, nudge_y  = 1, color="black") +
	# 	scale_color_distiller(palette = "Spectral")

	# MultiPhylo
	# dd_all = dd1 %>% mutate(.id = "a") %>% mutate(pvalue_extrinsic = as.numeric(pvalue_extrinsic)) %>% bind_rows(
	# 	dd2  %>% mutate(.id = "b")
	# ) %>% mutate(.id = as.factor(.id)) %>% mutate(estimate_extrinsic = rnorm(n()))
	#
	# trees <- list(tt_phylo1, tt_phylo2)
	# names(trees) = c("a", "b")
	# class(trees) <- "multiPhylo"
	# ggtree(trees)   + facet_wrap(~.id, scale="free") + geom_tiplab()
	# ggtree(trees, aes(color=estimate_extrinsic))  %<+% dd_all  + facet_wrap(~.id, scale="free") + geom_tiplab()
	#

	(
		switch(
			(!length(na.omit(dd$estimate_extrinsic))>1) + 1,

			tt_phylo %>% ggtree(	aes(color=estimate_extrinsic, size=`Cell type proportion`), layout="circular"	) %<+% dd +
			  scale_color_distiller(palette = "Spectral", na.value = 'gray87', direction = 1),

			tt_phylo %>% ggtree(	aes(color="Non significant", size=`Cell type proportion`), layout="circular"	) %<+% dd +
				scale_color_manual(values= c("Non significant" = "gray87" ))
		) +
		xlim(-2, NA) +
		geom_tiplab(size=4,  color="black", aes(angle=0),   offset=10, hjust = 0.5) +
		geom_nodelab(size=4,  color="black",  nudge_x = -20) +
		labs(color='Significant direction')
	) %>%
	rotate_tree(-154) +  theme(legend.position = "top" , legend.direction = "horizontal" )

	}


#' Create polar plot of results
#' @rdname ARMET_plotTree
#'
#' Prints a report of the hipothesis testing
#'
#' @import ggplot2
#' @import tibble
#' @import dplyr
#' @import data.tree
#' @importFrom ape as.phylo
#'
#'
#' @param ARMET-tc object
#'
#' @return a ggplot
#'
#' @export
ARMET_plotPolar = function(
	obj,
	size_geom_text = 3.5,
	my_breaks=c(0, 0.01, 0.03, 0.05, 0.1,0.3,0.5,0.7,1),
	prop_filter = 0.005,
	barwidth = 0.5,
	barheight = 2,
	legend_justification = 0.67,
	fill_direction = 1
){

	# Build phylo tree
	tt = obj$stats
	tt_phylo = tt %>% data.tree::as.phylo.Node()

	# Formatted names
	ct_names_formatted =
		ToDataFrameTree(tt, "name") %>%
		as_tibble %>% select(-levelName) %>%
		bind_cols(taxa = c(
			"Root",	"Epi", "Endo", "Fibro", "Immune", "Granulo", "Eosin", "Neutro", "Mono deriv", "Mono", "M0", "M1", "M2", "Dendr rest", "Drondr act", "Mast", "Activ", "Rest", "B", "Naive","Mem", "T", "CD8", "H1", "H2", "Follic", "γδ", "Reg", "M. cen", "M. act", "NK", "Activ", "Rest", "Plasma"
		))

	# Give formatted names
	tt_phylo$tip.label = ct_names_formatted %>% filter(name %in% tt_phylo$tip.label) %>% arrange(match(name, tt_phylo$tip.label)) %>% pull(taxa)
	tt_phylo$node.label = ct_names_formatted %>% filter(name %in% tt_phylo$node.label) %>% arrange(match(name, tt_phylo$node.label)) %>% pull(taxa)

	# Create annotation
	internal_branch_length = 40
	external_branch_length = 10
	dd =
		ToDataFrameTree(
			tt,
			"name",
			"isLeaf", "level", "leafCount",
			"Estimate",
			"Sig",
			"Cell type proportion",
			"Driver"
		) %>%
		as_tibble() %>%

		# BUG WITH EXCLUDED CELL TYOES
		filter(`Cell type proportion` %>% is.na %>% `!`) %>%

		left_join(ct_names_formatted, by= "name") %>%
		select(-levelName) %>% select(taxa, everything()) %>%
		mutate(Estimate = ifelse(Sig == "*" | Driver == "*", Estimate, NA)) %>%
		mutate(branch_length = ifelse(isLeaf, 0.1, 2)) %>%

		# Correct branch length
		mutate(level = level -1) %>%
		mutate(branch_length = ifelse(!isLeaf, internal_branch_length,	external_branch_length) ) %>%
		mutate(branch_length =  ifelse(	isLeaf, branch_length + ((max(level) - level) * internal_branch_length),	branch_length	))

	xx = dd %>%
		filter(taxa!="Adipo") %>%
	 	mutate(leafCount_norm = leafCount/max(leafCount)) %>%
		group_by(level) %>%
		arrange(desc(leafCount>1)) %>%
		mutate(leafCount_norm_cum = cumsum(leafCount_norm)) %>%
		ungroup() %>%
		arrange(level) %>%
		mutate(length_error_bar = leafCount_norm - 0.005) %>%
		mutate(x = leafCount_norm_cum - 0.5 * leafCount_norm) %>%
		mutate(angle = x * -360) %>%
		mutate(angle = ifelse(angle < - 240 | angle > -120, angle, angle - 180)) %>%
		mutate(angle = ifelse(angle > - 360, angle, angle + 360)) %>%
		mutate(`Cell type proportion` = ifelse(level<max(level), `Cell type proportion` + (0.008*level), `Cell type proportion`  )) %>%
		filter(level > 0) %>%
		mutate(`Cell type proportion norm` = `Cell type proportion` / max(`Cell type proportion`))

	my_rescale = function(x) { as.character( format(round( x * xx %>% pull(`Cell type proportion`) %>% max, 2), nsmall = 2)) }
	cusotm_root_trans = function() scales::trans_new("cusotm_root",function(x) x^(1/4), function(x) x^(4))
	DTC_scale =  xx %>% pull(Estimate) %>% abs() %>% { switch( (!all(is.na(.))) + 1, 1, .) } %>% max(na.rm = T)
	y_text_out = 5.0625
	y_text_in = 1.4641


	xx %>%
	{
		# Case if none is significant
		switch(
			(!length(na.omit(xx$Estimate))>0) + 1,

			# If there are significant
			ggplot(data=(.), aes(x = x,fill = Estimate,size = 1/sqrt(level))),

			# If there are not
			ggplot(data=(.),aes(x = x,fill = "Non significant",size = 1/sqrt(level)))
		)
	}	+
		annotate("rect", xmin=0, xmax=1, ymin=1, ymax= 2.0736, fill="grey95") +
		annotate("rect", xmin=0, xmax=1, ymin=2.0736, ymax=3.8416, fill="grey90") +
		annotate("rect", xmin=0, xmax=1, ymin=3.8416, ymax=6.5536, fill="grey87") +
	geom_bar(
		data = xx %>% filter(level == 1),
		aes(width = leafCount_norm, y = `Cell type proportion norm`), # / leafCount_norm),
		color = "grey20",  stat = "identity"
	) +
	geom_errorbar(
		data = xx %>% filter(level == 1) %>% mutate(x = ifelse(taxa=="Immune", 0.5, x)),
		aes(width = 0 , ymin=`Cell type proportion norm`, ymax=y_text_out-0.7), # / leafCount_norm),
		color = "grey20",  stat = "identity"
	) +
	geom_text(
		data =
			xx %>%
			filter(level == 1) %>%
			mutate(x = ifelse(taxa=="Immune", 0.5, x)) %>%
			mutate(angle = ifelse(taxa=="Immune", -0, angle)) %>%
			mutate(taxa = paste(taxa, Driver)) %>%
			mutate(taxa = sub("\\s+$", "", taxa)),
		aes(label=taxa, y = y_text_out, angle= angle  ) ,size =size_geom_text ) +
	#scale_x_continuous(labels = xx %>% filter(level == 1) %>% pull(taxa), breaks = xx %>% filter(level == 1) %>% pull(leafCount_norm_cum) - 0.5 * xx %>% filter(level == 1) %>% pull(leafCount_norm)) +
	geom_bar(
			data = xx %>% filter(level == 2),
			aes(width = leafCount_norm, y = `Cell type proportion norm` ), # / leafCount_norm),
			color = "grey20",  stat = "identity"
		) +
	geom_errorbar(
		data = xx %>% filter(level == 2),
		aes(width = 0 , ymin=`Cell type proportion norm`, ymax=2.2736), # / leafCount_norm),
		color = "grey20",  stat = "identity",linetype="dotted"
	) +
	geom_text(
		data =
			xx %>%
			filter(level == 2) %>%
			mutate(taxa = paste(taxa, Driver)) %>%
			mutate(taxa = sub("\\s+$", "", taxa)),
		aes(label=taxa, y = 2.8561, angle= angle) ,size =size_geom_text) +
	#scale_x_continuous(labels = xx %>% filter(level == 2) %>% pull(taxa), breaks = xx %>% filter(level == 2) %>% pull(leafCount_norm_cum) - 0.5 * xx %>% filter(level == 2) %>% pull(leafCount_norm)) +
	geom_bar(
		data = xx %>% filter(level == 3),
		aes(width = leafCount_norm, y = `Cell type proportion norm` ), # / leafCount_norm),
		color = "grey20", stat = "identity"
	)  +
	{
		# Make plotting robust if no level 3 cell types were detected
		switch(
			(!  xx %>% filter(level == 3) %>% filter(`Cell type proportion norm`>prop_filter | !is.na(Estimate)) %>% nrow() > 0) + 1,

			# If there are cell types
			geom_errorbar(
				data = xx %>% filter(level == 3) %>% filter(`Cell type proportion norm`>prop_filter | !is.na(Estimate)),
				aes(width = 0 , ymin=`Cell type proportion norm`, ymax=y_text_in-0.2), # / leafCount_norm),
				color = "grey20",  stat = "identity",linetype="dotted"
			) ,

			# If there are NOT cell types
			geom_errorbar(ymin=0, ymax=0)
		)
	} +
	{
		# Make plotting robust if no level 3 cell types were detected
		switch(
			(!  xx %>% filter(level == 3) %>% filter(`Cell type proportion norm`>prop_filter | !is.na(Estimate)) %>% nrow() > 0) + 1,

			# If there are cell types

				geom_text(
					data =
						xx %>%
						filter(level == 3)  %>%
						filter(
							`Cell type proportion norm`>prop_filter |
							!is.na(Estimate)
						) %>%
						mutate(taxa = paste(taxa, Driver)) %>%
						mutate(taxa = sub("\\s+$", "", taxa)),
					aes(label=taxa, y = y_text_in , angle= angle) ,size =size_geom_text
				),

			# If there are NOT cell types
			geom_errorbar(ymin=0, ymax=0)
		)
	} +
	{
		# Case if none is significant
		switch(
			(!length(na.omit(xx$Estimate))>0) + 1,

			# If there are significant
			scale_fill_distiller(
				palette = "Spectral",
				na.value = 'white',
				direction = fill_direction, name = "Trend",
				limits=c(-DTC_scale, DTC_scale)
			) ,

			# If there are not
			scale_fill_manual(values= c("Non significant" = "white" ), guide=FALSE)
		)

	} +
	{
		# Case if none is significant
		switch(
			(!length(na.omit(xx$Estimate))>0) + 1,

			# If there are significant
				guides(
					fill = guide_colorbar(
						label.position = "left",
						title.position = "left",
						title.hjust = 0.5,
						override.aes=list(fill=NA),
						ticks.colour = "black",
						barwidth = barwidth,
						barheight = barheight
					)
				) ,

			# If there are not
			scale_fill_manual(values= c("Non significant" = "white" ), guide=FALSE)
		)

	} +
	scale_y_continuous( breaks=my_breaks, trans = cusotm_root_trans(), labels = my_rescale ) +
	scale_size(range = c(0.3, 0.8), guide=FALSE) +
	theme_bw() +
	theme(
		axis.text.x = element_blank(),
		axis.ticks.x.top = element_blank(),
		axis.title.x=element_blank(),
		axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0), hjust = 0.65, vjust = 1),
		panel.border = element_blank(),
		axis.line.y = element_line(),
		panel.grid  = element_blank(),
		legend.position=c(0,0.05),
		legend.justification=c(legend_justification, 0),
		legend.title=element_text(angle = 90),
		legend.background = element_rect(colour = "transparent", fill = alpha("red", 0))
	)	+ ylab("Cell type proportion") +
		coord_polar(theta = "x")





}


ARMET_scatter = function(model.in, sample_idx=1){

	model.in$y %>%
		mutate(sample_idx = 1:n()) %>%
		gather(gene, observed, -sample_idx) %>%
		inner_join(
			model.in$x %>% bind_cols(model.in$x_genes),
			by="gene"
	) %>%
		gather(ct, reference, 4:7) %>%
		filter(sample_idx%in%!!sample_idx) %>%
		{
		ggplot(data=(.),aes(x=reference+1, y=observed+1, color=ct, gene=gene)) +
				geom_point() +
				scale_y_log10() +
				scale_x_log10()
		} %>% ggplotly()
}
