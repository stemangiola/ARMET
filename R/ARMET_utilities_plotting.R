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
