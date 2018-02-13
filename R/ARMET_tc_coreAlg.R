# ARMET_tc_coreAlg.R

ARMET_tc_coreAlg = function(
	obj.in, 
	ct, 
	real_prop_obj=NULL, 
	is_test=F){

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
	bg_tree =               obj.in$bg_tree
	sigma_hyper_sd =        obj.in$sigma_hyper_sd
	phi_hyper_sd =          obj.in$phi_hyper_sd
	alpha_hyper_value =     obj.in$alpha_hyper_value
	save_report =           obj.in$save_report

	# Calculate proportion of 0s
	theta = apply(mix, 2, function(mc) length(which(mc>0)))
	theta[theta==nrow(mix)] = theta[theta==nrow(mix)]-1
	theta = theta/nrow(mix) 
	theta = t(as.matrix(theta))
	
	ref.obj =   node_to_balanced_sampling(my_tree, ct, ref)
	ref =       cbind(do.call("cbind", lapply(ref.obj$ct_main, function(ro) ro$df)))
	markers =   unique(unlist(lapply(ref.obj$ct_main, function(ro) get_genes(node_from_name(my_tree, ro$name), recursive = F))))
	ref.cell_types = factor(colnames(ref))
	order_cell_types = levels(ref.cell_types)
	if(any(!markers%in%rownames(mix)))
		warning(sprintf("ARMET: The following markers were not present into the mixture dataset: %s Please reconsider the filtering process.\n", paste(markers[!markers%in%rownames(mix)])))
	
	markers =                           markers[markers%in%rownames(ref) & markers%in%rownames(mix)]
	mix =                               mix[markers,, drop=FALSE]
	ref =                               ref[markers,, drop=FALSE]
	ref.mean =                          get_mean_signature(ref)
	ref.mean =                          ref.mean[,order_cell_types, drop=FALSE]
	e.obj =                             prepare_input(ref, order_cell_types, markers)
	
	#/
	#|--------------------------------------------------------------------------
	#| Background
	#|--------------------------------------------------------------------------
	#|
	#| If down the my_tree compile the background\
	#\
	
	h = rev(rev(get_hierarchy(my_tree, ct))[-1])


	bg_trees = divide_trees_proportion_across_many_trees(bg_tree)
	bg_trees = lapply(bg_trees, add_absolute_proportions_to_tree)
	
	# 																										
	# lapply(names(bg_trees), function(n) get_proportions_array(bg_trees[n], n))
	# 
	# get_node_label_level_specfic(bg_trees[[1]], label = "absolute_proportion", recursive = T, level= 0, start_level = 0, stop_level = 0)
	# 
	if(length(h) > 0) {
		
		beta_bg = t(as.data.frame(lapply(bg_trees, get_last_existing_leaves_with_annotation)))
		bg = cbind(do.call("cbind", lapply(ref.obj$ct_background, function(ro) ro$df)))
		bg = bg[markers,,drop=F]
		
	} else {
		bg = matrix(rep(0, length(markers)))
		rownames(bg) = markers
		colnames(bg) = "nn"
		
		beta_bg = data.frame(rep(1, ncol(mix)), rep(0, ncol(mix))) 
		colnames(beta_bg) = c(ct, "bg")
		
	}
	
	p_target = t(beta_bg[,ct, drop=FALSE])
	beta_bg = as.matrix(beta_bg[,!colnames(beta_bg)%in%ct])
	
	# R bug for t() for 1 colmn matrices
	if(ncol(mix)==1) beta_bg = t(beta_bg)

	x_bg = get_mean_signature(bg)
	bg.obj = prepare_input(bg, levels(factor(colnames(bg))), markers)

	# Reorder the array same as the signatures
	if(ncol(beta_bg)>1) beta_bg = beta_bg[,colnames(x_bg),drop=F]

	beta_bg = as.matrix(beta_bg)
	
	# R bug for t() for 1 colmn matrices
	if(ncol(mix)==1) beta_bg = t(beta_bg)
	
	y_hat_background = as.matrix(beta_bg) %*% t(x_bg)
	rownames(y_hat_background) = colnames(mix)
	
	# balance the evidences
	e_map_matrix_dim = table(e.obj$e$ct_num)
	e_map_matrix_norm = median(e_map_matrix_dim)/e_map_matrix_dim
	e_i_matrix = do.call(plyr:::rbind.fill, lapply(unique(sort(e.obj$e$ct_num)), function(ct) data.frame(t(data.frame(rownames(e.obj$e[e.obj$e$ct_num==ct,] ))))))
	e_i_matrix = apply(e_i_matrix, 2, function(mc) {
		x = as.numeric(as.character(mc))
		x[is.na(x)] = 0
		x
	})
	
	# if design is NULL
	#my_design = NULL
	my_local_design = if(is.null(my_design)) matrix(rep(1, ncol(mix)), ncol=1) else my_design
	
	#/
	#|--------------------------------------------------------------------------
	#|--------------------------------------------------------------------------
	#\
	
	# Create input object for the model
	model.in = list(
		G =                               length(markers),
		S =                               ncol(mix),
		P =                               length(unique(e.obj$e$ct_num)),
		R =                               ncol(my_local_design),
		y =                               t(mix), # - mean(exp(t(mix)))) / sd(exp(t(mix))),
		X =                               my_local_design,
		x =                               ref.mean,
		y_hat_background =                y_hat_background,
		p_target =                        p_target,
		theta =                           theta,
		sigma_hyper_sd =                  sigma_hyper_sd,
		phi_hyper_sd =                    phi_hyper_sd,
		alpha_hyper_value =               alpha_hyper_value,
		
		
		is_mix_microarray =               as.numeric(is_mix_microarray),
		
		# For full Bayesian
		
		# Main node
		E =                               nrow(e.obj$e),
		E_MU =                            nrow(e.obj$e_mu),
		map_e_genes =                     e.obj$e$gene_num,
		map_e_ct =                        e.obj$e$ct_num,
		map_e_to_mu =                     e.obj$e$map_to_mu,
		e_ =                              e.obj$e$value,
		map_e_mu_gene =                   e.obj$e_mu$gene_num,
		map_e_mu_ct =                     e.obj$e_mu$ct_num,
		
		e_i_matrix =                      e_i_matrix,
		e_map_matrix_dim =                e_map_matrix_dim,
		e_map_matrix_norm =               e_map_matrix_norm,
		
		
		# Background
		B =                               nrow(bg.obj$e),
		B_MU =                            nrow(bg.obj$e_mu),
		Q =                               length(unique(bg.obj$e$ct_num)),
		map_bg_genes =                    bg.obj$e$gene_num,
		map_bg_ct =                       bg.obj$e$ct_num,
		map_bg_to_mu =                    bg.obj$e$map_to_mu,
		b =                               bg.obj$e$value,
		map_bg_mu_gene =                  bg.obj$e_mu$gene_num,
		map_bg_mu_ct =                    bg.obj$e_mu$ct_num,
		beta_bg =                         beta_bg
		
	)

	if(save_report) save(model.in, file=sprintf("%s/%s_model_in.RData", output_dir, ct))

	# Choose model
	model = if(fully_bayesian) stanmodels$ARMET_tc_recursive else stanmodels$ARMET_tcFix_recursive

	# Run model
	fit = 
		rstan::sampling(
			model,
			data=                             model.in,
			iter=                             1000 ,
			control =                         list(adapt_delta = 0.99, stepsize = 0.01, max_treedepth =15),
			cores=4
		)
	
	
	# Parse results

	proportions =                       t(parse_summary_vector_in_2D(apply( as.matrix(fit, pars = "beta"), 2, mean)))
	proportions_sd =                    t(parse_summary_vector_in_2D(apply( as.matrix(fit, pars = "beta"), 2, sd)))
	proportions_2.5 =                   t(parse_summary_vector_in_2D(apply( as.matrix(fit, pars = "beta"), 2, quantile, 0.025)))
	proportions_97.5 =                  t(parse_summary_vector_in_2D(apply( as.matrix(fit, pars = "beta"), 2, quantile, 0.975)))
	rownames(proportions) = rownames(proportions_2.5) = rownames(proportions_97.5) = rownames(proportions_sd) =  colnames(mix)
	colnames(proportions) =  colnames(proportions_2.5) = colnames(proportions_97.5) = colnames(proportions_sd) = order_cell_types

	# Save output
	if(save_report) save(fit, file=sprintf("%s/%s_fit.RData", output_dir, ct))
	p = rstan:::traceplot(fit, pars=c("alpha"), inc_warmup=F)
		
	#if(save_report) ggplot2:::ggsave(sprintf("%s_chains.png", ct), p)
	if(save_report)	write.csv(proportions, sprintf("%s/%s_composition.csv", output_dir, ct))
	
	
	proportions_4_plot = as.data.frame(proportions)
	proportions_4_plot$sample = rownames(proportions_4_plot)
	proportions_4_plot$cov = if(length(cov_to_test)>0) my_design[,"cov"] else 1
	proportions_4_plot = reshape:::melt(proportions_4_plot, id.vars=c("sample", "cov"))
	proportions_4_plot$perc_1 = reshape:::melt(proportions_2.5)$value
	proportions_4_plot$perc_2 = reshape:::melt(proportions_97.5)$value
	proportions_4_plot$minus_sd = proportions_4_plot$value - reshape:::melt(proportions_sd)$value
	proportions_4_plot$plus_sd = proportions_4_plot$value + reshape:::melt(proportions_sd)$value
	p = ggplot2:::ggplot( proportions_4_plot, ggplot2:::aes(x=jitter(cov), y=value,fill=factor(variable)))+ 
		ggplot2:::geom_boxplot(coef = 6) + 
		ggplot2:::geom_point(size=0.1) +
		ggplot2:::geom_linerange(ggplot2:::aes(ymin = perc_1, ymax = perc_2),  alpha=0.05) +
		ggplot2:::geom_errorbar(ggplot2:::aes(ymin = minus_sd, ymax = plus_sd),width=.2,  alpha=0.2) +
		ggplot2:::facet_grid(~variable) +
		ggplot2:::theme(
			axis.text.x= ggplot2:::element_text(angle=90, vjust=0.4,hjust=1), 
			panel.background = ggplot2:::element_blank()
		)
	
	#if(save_report) ggplot2:::ggsave(sprintf("%s/%s_proportions.pdf", output_dir, ct), useDingbats=FALSE,plot=p)
	
	if(length(cov_to_test)>0){
		#plot generate quantities
		gq <- as.matrix(fit, pars = c("beta_gen"))[1,]
		gq= parse_summary_vector_in_2D(gq)
		colnames(gq) = colnames(mix)
		rownames(gq) = colnames(model.in$x)
		gq = as.data.frame(t(gq))
		gq$cov = my_design[,"cov"]
		gq$sample = rownames(gq)
		gq = reshape:::melt(gq, id.vars=c("sample", "cov"))
		p = ggplot2:::ggplot( gq, ggplot2:::aes(x=factor(cov), y=value,fill=factor(variable)))+
			ggplot2:::geom_boxplot()+ 
			ggplot2:::geom_jitter(position=position_dodge(width=0.75), ggplot2:::aes(group=variable)) +
			ggplot2:::facet_grid(~variable) + 
			ggplot2:::theme(axis.text.x=ggplot2:::element_text(angle=90, vjust=0.4,hjust=1))
		
		#if(save_report) ggplot2:::ggsave(sprintf("%s/%s_proportions_posterior.pdf", output_dir, ct), useDingbats=FALSE,plot=p)
	}
	# Calculate predicted values
	#pred.df=as.matrix(ref[markers,])%*%t(as.matrix(res.df))
	
	
	# about means
	df_4_plot = data.frame(
		merge(
			reshape:::melt(mix+1), 
			reshape:::melt(t(proportions %*% t(ref.mean)+1)), 
			by=c("X1", "X2")
		)
	)
	colnames(df_4_plot)[3:4] = c("value", "predicted")
	df_4_plot$ct = sapply(df_4_plot$X1, function(x1) node_from_gene(my_tree, x1)$name )
	df_4_plot$bg = apply(df_4_plot, 1, function(mr) 	y_hat_background[mr[2], mr[1]] )
	df_4_plot$bg_prop = df_4_plot$bg /df_4_plot$value
	df_4_plot$bg_prop[df_4_plot$bg_prop >10] = 10
	
	
	df_4_plot = df_4_plot[df_4_plot$X2%in%unique(df_4_plot$X2)[1:20], ]
	#df_4_plot = df_4_plot[grep("CD8", df_4_plot$X2), ]
	
	p = ggplot2:::ggplot(df_4_plot, ggplot2:::aes(predicted, value, label = X1, size=bg_prop,color = factor(ct))) + 
		ggplot2:::geom_point(alpha=0.5) + 
		ggplot2:::scale_size(range = c(0, 5)) + 
		ggplot2:::expand_limits(x = 0.1, y = 0.1) +
		ggplot2:::geom_abline(intercept=0, slope=1) + 
		ggplot2:::scale_x_log10() + 
		ggplot2:::scale_y_log10() +
		ggrepel:::geom_text_repel(
			ggplot2:::aes(color = ct),
			size = 2,
			segment.alpha= 0.2) +
		ggplot2:::facet_wrap( ~ X2, nrow=5) +
		ggplot2:::coord_fixed() +
		ggplot2:::theme(
			panel.background = ggplot2:::element_blank(), 
			legend.text=ggplot2:::element_text(size=30)
		)
	
	#if(save_report) ggplot2:::ggsave(sprintf("%s/%s_predicted_values.pdf", output_dir, ct), useDingbats=FALSE, width = 80, height = 80, units = "cm", plot=p)
	
	
	if(!is.null(observed_prop) & 0) {
		writeLines("ARMET: start plotting expected gene vallues")
		link = tibble:::as_tibble(do.call("rbind", lapply(unique(colnames(observed_prop)), function(op) {
			ct = get_hierarchy(my_tree, ct=op)
			c(op, ct[ct%in%colnames(ref.mean)])
		})))
		if(ncol(link)>=2){
			colnames(link) = c("ct", "ct_new")
			observed_prop = tibble:::as_tibble(reshape:::melt(observed_prop))
			colnames(observed_prop) = c("sample", "ct", "value")
			
			observed_prop = 
				dplyr:::left_join(observed_prop, link, by="ct") %>% 
				dplyr:::select(-ct) %>% 
				dplyr:::group_by(sample, ct_new) %>% 
				dplyr:::summarise(value = sum(value, na.rm=TRUE)) %>% 
				tidyr:::spread(ct_new, value)
			
			for(r in colnames(ref.mean)[!colnames(ref.mean)%in%colnames(observed_prop)]) {
				my_df = matrix(rep(0, nrow(observed_prop)))
				colnames(my_df) = r
				observed_prop = observed_prop %>% dplyr:::bind_cols(as_tibble(my_df)) 
			}
			
			observed_prop = as.data.frame(observed_prop)
			rownames(observed_prop) = observed_prop[,1]
			observed_prop = observed_prop[,-1]
			observed_prop = as.matrix(observed_prop)
			observed_prop = observed_prop[,colnames(ref.mean)]
			
			print(proportions)
			
			df_4_plot = data.frame(merge(reshape:::melt(mix+1), reshape:::melt(t(observed_prop %*% t(ref.mean)+1)), by=c("X1", "X2")))
			colnames(df_4_plot)[3:4] = c("value", "predicted")
			df_4_plot$ct = sapply(df_4_plot$X1, function(x1) node_from_gene(my_tree, x1)$name )
			df_4_plot$bg = apply(df_4_plot, 1, function(mr) 	y_hat_background[mr[2], mr[1]] )
			df_4_plot$bg_prop = df_4_plot$bg /df_4_plot$value
			df_4_plot$bg_prop[df_4_plot$bg_prop >10] = 10
	
			p = ggplot2:::ggplot( df_4_plot[df_4_plot$X2%in%unique(df_4_plot$X2)[1:20], ], ggplot2:::aes(predicted, value, label = X1, size=bg_prop)) + 
				ggplot2:::geom_point(alpha=0.5) + 
				ggplot2:::scale_size(range = c(0, 5)) + 
				ggplot2:::expand_limits(x = 0.1, y = 0.1) +
				ggplot2:::geom_abline(intercept=0, slope=1) + 
				ggplot2:::scale_x_log10() + 
				ggplot2:::scale_y_log10() +
				ggrepel:::geom_text_repel(
					aes(color = ct),
					size = 2,
					segment.alpha= 0.2) +
				ggplot2:::facet_wrap( ~ X2, nrow=5) +
				ggplot2:::coord_fixed() +
				ggplot2:::theme(
					panel.background = ggplot2:::element_blank(), 
					legend.text=ggplot2:::element_text(size=30)
				)
			
			#if(save_report) ggplot2:::ggsave(sprintf("%s/%s_predicted_values_observed_prop.pdf", output_dir, ct), useDingbats=FALSE, width = 80, height = 80, units = "cm", plot=p)
		}
		
	}
	
	
	list(proportions = proportions)
}
