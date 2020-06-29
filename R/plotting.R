#' Create polar plot of results
#' @rdname plot_polar
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
plot_polar = function(	.data,
											 size_geom_text = 3.5,
											 my_breaks=c(0, 0.01, 0.03, 0.05, 0.1, 0.3, 0.5, 0.7, 1),
											 prop_filter = 0.005,
											 barwidth = 0.5,
											 barheight = 2,
											 legend_justification = 0.67,
											 fill_direction = 1){
	
	xx  = 
		.data %>%
		calculate_x_for_polar %>%
		
		# Add angle text
		mutate(angle = x * -360) %>%
		mutate(angle = ifelse(angle < - 240 | angle > -120, angle, angle - 180)) %>%
		mutate(angle = ifelse(angle > - 360, angle, angle + 360)) %>%
		#mutate(median_proportion = ifelse(level<max(level), median_proportion + (0.008*level), median_proportion  )) %>%
		filter(level > 0) %>%
		
		# Calculate median proportions
		filter(.value_alpha1 %>% is.na %>% `!`) %>%
		mutate(summarised_proportion =
					 	proportions %>% map(~ .x %>% summarise(median_proportion = .value %>% median))
		) %>%
		unnest(summarised_proportion) %>%
		#mutate(median_proportion = median_proportion / max(median_proportion)) %>%
		
		distinct() %>%
		left_join(
			ct_names_polar, by=c("Cell type category" = "name")
		)
	
	my_rescale = function(x) { as.character( format(round( x * xx %>% pull(median_proportion) %>% max, 2), nsmall = 2)) }
	cusotm_root_trans = function() scales::trans_new("cusotm_root",function(x) x^(1/4), function(x) x^(4))
	DTC_scale =  xx %>% pull(Estimate) %>% abs() %>% { switch( (!all(is.na(.))) + 1, 1, .) } %>% max(na.rm = T)
	y_text_out = 8.5
	y_text_in = 1.4641
	
	# Size outer circles
	soc = (c(0, 0.2, 0.4, 0.6, 0.8)*0.5 + 1)^4  
	
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
		annotate("rect", xmin=0, xmax=1, ymin=soc[1], ymax= soc[2], fill="grey95") +
		annotate("rect", xmin=0, xmax=1, ymin=soc[2], ymax=soc[3], fill="grey90") +
		annotate("rect", xmin=0, xmax=1, ymin=soc[3], ymax=soc[4], fill="grey87") +
		annotate("rect", xmin=0, xmax=1, ymin=soc[4], ymax=soc[5], fill="grey83") +
		
		# Lv 1
		geom_bar(
			data = xx %>% filter(level == 1),
			aes(width = leafCount_norm, y = `median_proportion`), # / leafCount_norm),
			color = "grey20",  stat = "identity"
		) +
		geom_errorbar(
			data = xx %>% filter(level == 1), # %>% mutate(x = ifelse(`Cell type category`=="immune_cell", 0.5, x)),
			aes(width = 0 , ymin=`median_proportion`, ymax=mean(soc[1:2])), # / leafCount_norm),
			color = "grey20",  stat = "identity"
		) +
		geom_text(
			data =
				xx %>%
				filter(level == 1) %>%
				# mutate(x = ifelse(`Cell type category`=="immune_cell", 0.5, x)) %>%
				# mutate(angle = ifelse(`Cell type category`=="immune_cell", -0, angle)) %>%
				mutate(`formatted` = sub("\\s+$", "", `formatted`)),
			aes(label=`formatted`, y = mean(soc[1:2]), angle= angle  ) ,size =size_geom_text ) +
		#scale_x_continuous(labels = xx %>% filter(level == 1) %>% pull(`formatted`), breaks = xx %>% filter(level == 1) %>% pull(leafCount_norm_cum) - 0.5 * xx %>% filter(level == 1) %>% pull(leafCount_norm)) +
		
		# Lv 2
		geom_bar(
			data = xx %>% filter(level == 2),
			aes(width = leafCount_norm, y = `median_proportion` ), # / leafCount_norm),
			color = "grey20",  stat = "identity"
		) +
		geom_errorbar(
			data = xx %>% filter(level == 2),
			aes(width = 0 , ymin=`median_proportion`, ymax=mean(soc[2:3])), # / leafCount_norm),
			color = "grey20",  stat = "identity",linetype="dotted"
		) +
		geom_text(
			data =
				xx %>%
				filter(level == 2) %>%
				mutate(`formatted` = sub("\\s+$", "", `formatted`)),
			aes(label=`formatted`, y = mean(soc[2:3]), angle= angle) ,size =size_geom_text) +
		#scale_x_continuous(labels = xx %>% filter(level == 2) %>% pull(`formatted`), breaks = xx %>% filter(level == 2) %>% pull(leafCount_norm_cum) - 0.5 * xx %>% filter(level == 2) %>% pull(leafCount_norm)) +
		
		# Lv 3
		geom_bar(
			data = xx %>% filter(level == 3),
			aes(width = leafCount_norm, y = `median_proportion` ), # / leafCount_norm),
			color = "grey20", stat = "identity"
		)  +
		{
			# Make plotting robust if no level 3 cell types were detected
			switch(
				(!  xx %>% filter(level == 3) %>% filter(`median_proportion`>prop_filter | !is.na(Estimate)) %>% nrow() > 0) + 1,
				
				# If there are cell types
				geom_errorbar(
					data = xx %>% filter(level == 3) %>% filter(`median_proportion`>prop_filter | !is.na(Estimate)),
					aes(width = 0 , ymin=`median_proportion`, ymax=mean(soc[3:4])), # / leafCount_norm),
					color = "grey20",  stat = "identity",linetype="dotted"
				) ,
				
				# If there are NOT cell types
				geom_errorbar(ymin=0, ymax=0)
			)
		} +
		{
			# Make plotting robust if no level 3 cell types were detected
			switch(
				(!  xx %>% filter(level == 3) %>% filter(`median_proportion`>prop_filter | !is.na(Estimate)) %>% nrow() > 0) + 1,
				
				# If there are cell types
				
				geom_text(
					data =
						xx %>%
						filter(level == 3)  %>%
						filter(
							`median_proportion`>prop_filter |
								!is.na(Estimate)
						) %>%
						mutate(`formatted` = sub("\\s+$", "", `formatted`)),
					aes(label=`formatted`, y = mean(soc[3:4]) , angle= angle) ,size =size_geom_text
				),
				
				# If there are NOT cell types
				geom_errorbar(ymin=0, ymax=0)
			)
		} +
		
		
		# Lv 4
		geom_bar(
			data = xx %>% filter(level == 4),
			aes(width = leafCount_norm, y = `median_proportion` ), # / leafCount_norm),
			color = "grey20", stat = "identity"
		)  +
		{
			# Make plotting robust if no level 4 cell types were detected
			switch(
				(!  xx %>% filter(level == 4) %>% filter(`median_proportion`>prop_filter | !is.na(Estimate)) %>% nrow() > 0) + 1,
				
				# If there are cell types
				geom_errorbar(
					data = xx %>% filter(level == 4) %>% filter(`median_proportion`>prop_filter | !is.na(Estimate)),
					aes(width = 0 , ymin=`median_proportion`, ymax=mean(soc[4:5])), # / leafCount_norm),
					color = "grey20",  stat = "identity",linetype="dotted"
				) ,
				
				# If there are NOT cell types
				geom_errorbar(ymin=0, ymax=0)
			)
		} +
		{
			# Make plotting robust if no level 4 cell types were detected
			switch(
				(!  xx %>% filter(level == 4) %>% filter(`median_proportion`>prop_filter | !is.na(Estimate)) %>% nrow() > 0) + 1,
				
				# If there are cell types
				
				geom_text(
					data =
						xx %>%
						filter(level == 4)  %>%
						filter(
							`median_proportion`>prop_filter |
								!is.na(Estimate)
						) %>%
						mutate(`formatted` = sub("\\s+$", "", `formatted`)),
					aes(label=`formatted`, y = mean(soc[4:5]) , angle= angle) ,size =size_geom_text
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

#' 
#' #' Create polar plot of results
#' @rdname plot_scatter
#'
#' Prints a report of the hipothesis testing
#'
#' @import ggplot2
#' @import tibble
#' @import dplyr
#'
#'
#' @param ARMET-tc object
#'
#' @return a ggplot
#'
#' @export
#' 
plot_scatter = function(.data){
	
	# data_CI = 
	# 	.data$proportions %>%
	# 	select(-proportions, -.variable) %>%
	# 	#unnest(draws) %>%
	# 	select( -draws) %>%
	# 	unnest(rng)  %>% 
	# 	group_by( C, Q, .variable, level, `Cell type category`, Rhat) %>% 
	# 	tidybayes::median_qi() %>%
	# 	ungroup() %>%
	# 	left_join(
	# 		.data$proportions %>%
	# 			select(-draws, -.variable) %>%
	# 			unnest(proportions) %>% distinct(Q,  DFS_MONTHS)
	# 	) %>%
	# 	select(level,  `Cell type category`, Q, .upper, .lower)
	
	inferred_y = 	
		.data$proportions %>%
		select(level, `Cell type category`, proportions) %>%
		unnest(proportions) %>% 
		select(level, `Cell type category`, Q, .draws) %>%
		unnest(.draws) %>%
		rename(inferred_y = .value_relative)
	
	
	inferred_x = 
		map2_dfr(
			.data$internals$fit,
			1:length(.data$internals$fit),
			~ .x %>%
				tidybayes::gather_draws(X_[Q, A]) %>% 
				#tidybayes::median_qi() %>% 
				ungroup() %>% 
				filter(A==2) %>% 
				mutate(level = .y)
		) %>%
		rename(inferred_x = .value) 
	
	plot_data = 
		inferred_y %>%
		left_join(inferred_x) %>%
		sample_frac(0.1) %>%
		left_join( .data$proportions %>% select(proportions) %>% unnest(proportions) %>% distinct(sample, Q, alive) ) %>%
		#filter(`Cell type category` == "epithelial") %>%
		#mutate(inferred_x = inferred_x + rnorm(n(), 0, 0.01) %>% abs) %>%
		filter(inferred_y %>% is.na %>% `!`) %>%
		group_by(`Cell type category`, sample, Q, level, alive) %>%
		summarise(
			xmin= quantile(inferred_x, 0.05), xmax = quantile(inferred_x, 0.95), x = mean(inferred_x),
			ymin= quantile(inferred_y, 0.05), ymax = quantile(inferred_y, 0.95), y = mean(inferred_y)
		) %>%
		ungroup() %>%
		mutate(area = (xmax-xmin) * (ymax-ymin))
	
	outlier_df =
		.data$proportions %>%
		filter(map_lgl(rng, ~!is.null(.)) ) %>% 
		select(level, `Cell type category`, C,rng) %>% 
		unnest(rng) %>% 
		group_by(level, `Cell type category`, C,Q, .variable) %>% 
		tidybayes::median_qi(.width = 0.95) %>%
		ungroup() %>%
		rename( .lower_rng = .lower, .upper_rng = .upper, .median_rng = .value)
	

	plot_data = 
		plot_data %>%
		left_join(outlier_df) %>%
		rowwise() %>%
		mutate(outlier = y %>% between(.lower_rng, .upper_rng) %>% `!`) %>%
		ungroup()
	
	ggplot(plot_data, aes(x, y)) +
		
		geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax =ymax, alpha = -area), data =  plot_data %>% filter(alive)) +
		geom_errorbar( aes(ymin =ymin, ymax = ymax, color = outlier), data = plot_data %>% filter(!alive)) +
		geom_point( data =  plot_data %>% filter(alive), color="blue", shape=".") +
		geom_line(aes(x, .upper_rng)) +
		geom_line(aes(x, .lower_rng)) +
		#geom_density_2d(bins=3, fill = after_stat(density), geom = "polygon") +
		#geom_point(aes(color = alive)) +
		facet_wrap(~`Cell type category`, scale="free") +
		scale_alpha(range = c(0.1, 0.3)) +
		
		theme_bw() +
		theme(
			panel.border = element_blank(),
			axis.line = element_line(),
			panel.grid.major = element_line(size = 0.2),
			panel.grid.minor = element_line(size = 0.1),
			text = element_text(size = 12),
			legend.position = "bottom",
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
	
}

#' This is a generalisation of ifelse that acceots an object and return an objects
#'
#' @import ggplot2
#'
#' @export
plot_markers  = function(result, level, S = NULL, cores = 20){
	
	library(multidplyr)
	
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

plot_heatmap = function(.data){
	
	.data %>%
		get_signatures %>%
		filter(A==2, .variable %in% c("alpha_1", "alpha_a")) %>%
		lower_triangular %>%
		separate(`Cell type category_1`, c("l1"), sep="_", extra = "drop", remove = F) %>%
		separate(`Cell type category_2`, c("l2"), sep="_", extra = "drop", remove = F) %>%
		
		mutate(
			label = paste(
				substr(`l2`, start = 0, stop = 3),
				substr(`l1`, start = 0, stop = 3),
				sep = "\n"
			)
		) %>%
		ggplot(aes( `Cell type category_1`, `Cell type category_2`, fill=prob, label=label)) +
		geom_tile() +
		geom_text(angle = 45) +
		#facet_wrap(~.variable, scale="free") +
		scale_fill_distiller(palette = "Spectral", na.value="transparent", limits = c(-1,1) ) +
		coord_fixed() +
		scale_x_discrete(position = "top") +
		theme(
			plot.background = element_rect(fill = "transparent",colour = NA),
			panel.grid=element_blank(),
			panel.background=element_blank(),
			panel.border = element_blank(),
			plot.margin = unit(c(0, 0, 0, 0), "npc"),
			axis.title.x = element_blank(),
			axis.title.y = element_blank(),
			axis.text = element_blank()
			#,
			#axis.text.x = element_text(angle = 90, hjust = 0)
		)
	
	# plot_2 <- ARMET_TCGA_result_hierarchical_beta %>%
	# 	get_signatures %>%
	# 	filter(A==2) %>%
	# 	nest(data = -.variable) %>%
	# 	mutate(lower_tri = map(data, ~.x %>% lower_triangular)) %>%
	# 	mutate(p = map(
	# 		lower_tri,
	# 		~.x %>%
	# 			ggplot(aes( `Cell type category_1`, `Cell type category_2`, fill=prob)) +
	# 			geom_tile() +
	# 			#facet_wrap(~.variable, scale="free") +
	# 			scale_fill_distiller(palette = "Spectral", na.value="transparent", limits = c(-1,1) ) +
	# 			coord_fixed() +
	# 			scale_x_discrete(position = "top") +
	# 			theme(
	# 				plot.background = element_rect(fill = "transparent",colour = NA),
	# 				panel.grid=element_blank(),
	# 				panel.background=element_blank(),
	# 				panel.border = element_blank(),
	# 				plot.margin = unit(c(0, 0, 0, 0), "npc"),
	# 				axis.title.x = element_blank(),
	# 				axis.title.y = element_blank(),
	# 				axis.text.x = element_text(angle = 90, hjust = 0)
	# 			)
	# 	))
	
}

