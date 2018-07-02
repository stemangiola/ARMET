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
			"Root",	"Epi", "Endo", "Fibro", "Adipo", "Immune", "Granulo", "Eosin", "Neutro", "Mono deriv", "Mono", "M0 macro", "M1 macro", "M2 macro", "Rest", "Activ", "Mast", "Activ", "Rest", "B", "Naive","Mem", "T", "CD8", "CD4 naive", "H1", "H2", "Follic", "γδ", "Reg", "Mem cent", "Mem activ", "NK", "Activ", "Rest", "Plasma"
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

	# Plot circular tree
	(
		tt_phylo %>% 
		ggtree(aes(color=estimate_extrinsic, size=`Cell type proportion`), layout="circular") %<+% dd +  
		xlim(-2, NA) +
		#geom_text(aes(label=label, angle=angle), vjust =-.2, na.rm = T) +
		#geom_tiplab(size=5,  color="black" , offset = 2 ) + 
		geom_tiplab(size=4,  color="black", aes(angle=0),   offset=10, hjust = 0.5) +
		#geom_tiplab2(size=4,  color="black", aes(angle=angle),  align=TRUE,  offset=8) +
		geom_nodelab(size=4,  color="black",  nudge_x = -20) 
		) %>%
	rotate_tree(-154) +  
		theme(legend.position = "top" , legend.direction = "horizontal" ) +
		scale_color_distiller(palette = "Spectral", na.value = 'gray87', direction = 1)
	#scale_color_continuous(na.value = 'grey', high='red', low='blue') 
	}
