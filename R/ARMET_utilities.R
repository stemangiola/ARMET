
#' Check the input for anomalies
check_input = function(mix, is_mix_microarray, my_design, cov_to_test, prior_sd, custom_ref, tree){
	# Check if duplicated names
	if(any(table(colnames(mix))>1)) stop("ARMET: you have duplicated column names in the mix data frame")
	
	# Check if microarray data
	if(is_mix_microarray) writeLines("ARMET: The input matrix is from microarray.")
	
	# Check if on Windows
	if(!.Platform$OS.type == "unix") stop("ARMET: currenly ARMET works for GNU/Linux systems only.")
	
	# Check if NA in mix
	if(any(is.na(as.vector(mix)))) stop("ARMET: NAs found in the query matrix")
	
	# This is how many conditions are in the study (e.g., treatment-vs-non-treatment)
	if(!is.null(my_design)) {
		writeLines("The design matrix is :")
		print(head(my_design))
	}
	
	# Check if mix have colnames and rownames
	if(is.null(colnames(mix)) | is.null(rownames(mix))) stop("ARMET: query data frame (mix expression) has to have colnames (sample names) and rownames (gene names)")
	
	# Set up the covariate to test if any
	# if(
	#   (is.null(cov_to_test) & ncol(my_design)>1) |
	#   !cov_to_test%in%colnames(my_design)
	# ) stop(sprintf("ARMET: you have to specify one or more covariate(s) to test for your design matrix among these: %s", paste(colnames(my_design), collapse=" ")))
	if(is.null(my_design) & !is.null(cov_to_test)) stop("ARMET: you have specified a covariate to test but no design matrix")
	
	# Check prior sd
	if(prior_sd<=0) stop("ARMET: prior_sd must be a positive number.")
	
	# Check custom reference
	if(!is.null(custom_ref) & !all(get_leave_names(tree, last_level = 1)%in%colnames(custom_ref))) stop("ARMET: Some cell types within the tree are absent in your custom reference provided") 
	
}

#' Get parameters from ini file if it exists
get_ini = function(){
	if(file.exists("ARMET.ini")) {
		pars =  ini:::read.ini("ARMET.ini")$ARMET
		writeLines("Importing parameters from .ini file..")
		for(n in names(pars)){
			if(pars[n]%in%c("T", "F", "TRUE", "FALSE"))
				pars[n] = as.logical(pars[n])
			assign(n, unlist(pars[n]))
		}
	}
}
#' Creates the directory needed for archive
#'
#' @param name A char
create_temp_result_directory = function(name="ARMET_results"){
	time_id = gsub(" |:|-", "_", Sys.time())
	output_dir = sprintf("%s_%s", name, time_id)
	dir.create(output_dir, showWarnings = FALSE)
	output_dir
}

#' Remove gene redundancy averagin duplicated genes
#'
#' @param ref A matrix
#' @param symbol A char vector
#' @return Non redundant matrix
average_duplicated_genes = function(ref, symbol){
	temp = data.frame(symbol = symbol, ref)
	temp = na.omit(temp)
	us = unique(symbol)
	
	cl = parallel:::makeCluster(20)
	parallel:::clusterExport(cl, "temp", envir = environment())
	parallel:::clusterEvalQ(cl, library("matrixStats"))
	
	ref.symbol = do.call("rbind", parallel:::parLapply(cl, us, function(tt) {
		matrixStats:::colMedians(as.matrix(temp[temp$symbol==tt,-1]))
	}))
	parallel:::stopCluster(cl)
	
	colnames(ref.symbol) = colnames(ref)
	rownames(ref.symbol) = us
	
	ref.symbol
}

#' Plots the densities distributions columnwise of a matrix
#'
#' @param df A matrix
#' @param color A char
#' @param fill A array of numbers
#' @param alpha A number
#' @param do_log A bool
#' @return A ggplot
plot_densities = function(df, color="0",  fill = "0", alpha = 0.30, do_log = T){
	
	if(do_log) df = df %>% dplyr:::mutate(value=value+0.1)
	p = 
		ggplot2:::ggplot(	df, ggplot2:::aes(value, group=sample, color = factor(color))) +
		ggplot2:::geom_line(stat="density", alpha=alpha) +
		ggplot2:::expand_limits(x=0.1) +
		ggplot2:::theme_bw() +
		ggplot2:::theme(
			panel.border = ggplot2:::element_blank(),
			panel.grid.major = ggplot2:::element_blank(),
			panel.grid.minor = ggplot2:::element_blank(),
			axis.line = ggplot2:::element_line(colour = "black")
		)
	
	if(do_log) p = p + ggplot2:::scale_x_log10()
	
	return( p )
}

#' Plots the densities from a melted data frame
#'
#' @param df_melted A matrix
#' @return A ggplot
plot_densities_from_melted = function(df_melted){
	
	p = ggplot2:::ggplot(df_melted, ggplot2:::aes(value+0.1, group=X2, color=source)) +
		ggplot2:::geom_line(stat="density", alpha=0.15) +
		ggplot2:::scale_x_log10() +
		ggplot2:::expand_limits(x=0.1) +
		ggplot2:::theme_bw() +
		ggplot2:::theme(
			panel.border = ggplot2:::element_blank(),
			panel.grid.major = ggplot2:::element_blank(),
			panel.grid.minor = ggplot2:::element_blank(),
			axis.line = ggplot2:::element_line(colour = "black")
		)
	return( p )
}

#' Plots the densities of ref and mix
#'
#' @param reference A matrix
#' @param obj A matrix
#' @return A ggplot
plot_densities_double = function(reference,obj){
	
	reference = as.matrix(reference)
	obj = as.matrix(obj)
	
	# Build data fraes for plotting
	reference.4plot = reference
	colnames(reference.4plot) = paste0("r", 1:dim(reference.4plot)[2])
	reference.4plot = data.frame(reshape:::melt(reference.4plot), source="reference")
	obj.4plot = obj
	colnames(obj.4plot) = paste0("r", 1:dim(obj.4plot)[2])
	obj.4plot = data.frame(reshape:::melt(obj.4plot), source="obj")
	
	# Build plot
	p = plot_densities_from_melted(rbind(reference.4plot, obj.4plot))
	
	return( p )
}


#' Plots the chains for checking convergence
#'
#' @param fit A stan object
#' @param do_show A bool
#' @return A ggplot
sanity_check_p = function(fit, do_show = F){
	
	plot_chains = rstan:::traceplot(fit, inc_warmup = FALSE, pars="alpha")
	if(do_show) plot(plot_chains)
	
	return( plot_chains )
}

error_if_log_transformed = function(x){
	if(length(x$value)>0) if(max(x$value, na.rm=T)<50) 
		stop("ARMET: The input was log transformed in: check_if_sd_zero_and_correct")
}

#' Checks if the standard deviation of a gene is 0
#'
#' @param df A matrix
check_if_sd_zero_and_correct = function(df, node){
	
	# Sanity check
	error_if_log_transformed(df)
	
	# Function omit zeros
	zero.omit = function(x) x[x>0]
	
	# Summarize counts
	df.summary = df %>%
		dplyr:::group_by(gene, ct) %>%
		dplyr:::summarise(log_sigma = sd(log(value+1)), log_avg = mean(log(value+1))) %>%
		dplyr:::ungroup() %>%
		dplyr:::mutate(to_recalculate = log_sigma==0) 
	
	find_closest_sd_uo_the_tree = function(my_ct, my_gene, tb){
		
		# Iterate upwards thoward root to get the sd from the closest group
		hierarchy = get_hierarchy(node, my_ct)
		ancestors = rev(hierarchy[1:(length(hierarchy)-1)])
		
		for(h in ancestors){
			
			# Get first descendant of ancestor including the ancestor 
			ct_to_consider = unlist(c(h, get_leave_label(node_from_name(node, h), recursive = F)))
			
			mean_log_sigma = tb %>% 
				dplyr:::filter(ct %in% ct_to_consider & gene == my_gene) %>%
				dplyr:::summarise(mean_log_sigma = mean(zero.omit(log_sigma)))
			
			if(mean_log_sigma>0) break
			
		}
		
		as.numeric(mean_log_sigma)
	}
	
	# Add closest sd
	df.summary = df.summary %>%
		dplyr:::rowwise() %>%
		dplyr:::mutate(
			log_sigma = 
				ifelse(
					log_sigma==0,
					find_closest_sd_uo_the_tree(ct, gene, df.summary),
					log_sigma
				)
		) %>%
		dplyr:::ungroup()
	
	if(nrow(df.summary %>% filter(to_recalculate))>0){
		writeLines("ARMET: One of your markers has 0 variance, this is not compatible with MCMC inference.")
		writeLines("The following genes will acquire a background variance.")
		print( df.summary %>% filter(to_recalculate) )
	}
	
	# Sample for sd == 0
	df %>%
		left_join(df.summary, by=c("gene", "ct")) %>%
		dplyr:::rowwise() %>%
		dplyr:::mutate(value = ifelse(to_recalculate, exp(rnorm(1, log_avg, log_sigma)), value)) %>%
		dplyr:::ungroup() %>%
		dplyr:::select(-log_sigma, -log_avg, -to_recalculate)

}

#' Calculate the norm factor with calcNormFactor from limma
#'
#' @param df A matrix
#' @param reference A reference matrix
#' @param cpm_theshold A number
#' @param prop A number
#' @return A list including the filtered data frame and the normalization factors
rnaseq_norm.calcNormFactor = function(df, reference = NULL, cpm_theshold = 0.5, prop = 3/4){
	
	if(max(df$value, na.rm = T) < 50) stop("ARMET: Both mixture and signatures have to be in count form, log tranformation detected")
	
	df.filt = df %>% 
	{ if("ct"%in%names(.)) dplyr:::select(-ct)	else .}	%>%
		tidyr:::spread(sample, value) %>% 
		tidyr:::drop_na() %>%
		dplyr:::do(
			(.) %>% 
				dplyr:::filter(
					rowSums(
						edgeR:::cpm(
							(.) %>% 
								dplyr:::select(-gene)
						) > cpm_theshold
					) >= 
						ceiling(ncol( 
							(.) %>% dplyr:::select(-gene)
						) * prop)
				)
		)
	
	list(
		nf = tibble:::tibble(
			sample = factor(colnames(df.filt %>% dplyr:::select(-gene))),
			nf = edgeR:::calcNormFactors(df.filt %>% dplyr:::select(-gene), refColumn=reference)
		),
		df = df.filt %>% 
			tidyr:::gather(sample, value, 2:ncol(.)) %>% 
			dplyr:::mutate(sample = factor(sample)) %>%
			{ 
				if("ct"%in%names(.)) 
					dplyr:::left_join(
						df %>% dplyr:::distinct(sample, ct), 
						by="sample"
					) 
				else 
					.
			}
	)
}

#' Normalize a RNA seq data set using rnaseq_norm.calcNormFactor
#'
#' @param df A matrix
#' @param reference A reference matrix
#' @param cpm_theshold A number
#' @param prop A number
#' @return A list including the filtered data frame and the normalization factors
rnaseq_norm = function(df, reference = NULL, cpm_theshold = 0.5, prop = 3/4){
	writeLines("Normalizing RNA-seq data with TMM")
	
	nf.obj = rnaseq_norm.calcNormFactor(df, reference, cpm_theshold, prop)
	nf = nf.obj$nf
	df.filt = nf.obj$df
	
	df %>% 
		dplyr:::group_by(sample) %>% 
		dplyr:::mutate(tot = sum(value, na.rm = T)) %>%
		dplyr:::ungroup() %>%
		dplyr:::left_join(nf, by="sample") %>%
		dplyr:::left_join(
			df.filt %>%
				dplyr:::group_by(sample) %>%
				dplyr:::summarise(tot_filt = sum(value, na.rm = T)),
			"sample"
		) %>%
		dplyr:::mutate(value = value / (tot_filt * nf) * max(tot)) %>%
		dplyr:::select(-tot, -nf, -tot_filt)
	
}

#' Normalize ref to match mix using TMM
#'
#' @param ref A matrix
#' @param mix A matrix
#' @return A list including the ref and mix normalized and a ggplot
rnaseq_norm_ref_mix = function(ref, mix){
	
	if(max(ref$value, na.rm = T) < 50 | max(mix$value, na.rm = T) < 50) stop("ARMET: Both objects have to be in real scale, log tranformation detected")
	
	# Normalize the mix
	mix = rnaseq_norm(mix)
	
	# Calclate normalization factors for mix and ref => two numbers
	nf.obj = rnaseq_norm.calcNormFactor(
		rbind(
			mix %>% 
				dplyr:::group_by(gene) %>% 
				dplyr:::summarise(value = median(value, na.rm=T)) %>% 
				dplyr:::mutate(sample="mix"),
			ref %>% 
				dplyr:::group_by(gene) %>% 
				dplyr:::summarise(value = median(value, na.rm=T)) %>% 
				dplyr:::mutate(sample="ref")
		), 
		1)
	
	nf = nf.obj$nf
	my_df = nf.obj$df
	
	tot_reference = (my_df %>%
									 	dplyr:::filter(sample == "mix") %>% 
									 	dplyr:::summarise(tot = sum(value, na.rm=T)))$tot 
	
	tot_other = (my_df %>%
							 	dplyr:::filter(sample == "ref") %>% 
							 	dplyr:::summarise(tot = sum(value, na.rm=T))
	)$tot 
	
	mix = mix %>%
		dplyr:::mutate(tot_ref = tot_reference) %>%
		dplyr:::mutate(nf = 
									 	(nf %>%
									 	 	dplyr:::filter(sample=="mix"))$nf
		) %>%
		dplyr:::mutate(tot = tot_ref) %>%
		dplyr:::mutate(value = value / (tot * nf) * tot_ref) %>%
		dplyr:::select(-tot, -tot_ref, -nf)
	
	ref = ref %>%
		dplyr:::mutate(tot_ref = tot_reference) %>%
		dplyr:::mutate(nf = 
									 	(nf %>%
									 	 	dplyr:::filter(sample=="ref"))$nf
		) %>%
		dplyr:::mutate(tot = tot_other) %>%
		dplyr:::mutate(value = value / (tot * nf) * tot_ref) %>%
		dplyr:::select(-tot, -tot_ref, -nf)
	
	p = plot_densities( 
		rbind( 
			ref %>% 
				dplyr:::mutate(color=1) %>%
				dplyr:::select(-ct),
			mix %>% 
				dplyr:::mutate(color=1)
		)
	)
	
	return(list(ref = ref, mix = mix,  plot = p))
}

#' Normalize array data of ref and mix
#'
#' @param ref A matrix
#' @param mix A matrix
#' @return A list including the ref and mix normalized and a ggplot
array_norm = function(ref, mix){
	
	writeLines("Normalizing Array data")
	log_ex = log(cbind(ref, mix)+1)
	log_ex.norm = limma:::normalizeBetweenArrays(log_ex)
	ref = exp(log_ex.norm[,1:dim(ref)[2]])
	mix = exp(log_ex.norm[,(dim(ref)[2]+1):dim(log_ex.norm)[2]])
	
	# Build plot
	p = plot_densities(cbind(ref,mix))
	
	return(list(ref=ref, mix=mix, plot=p))
}

#' Normalize array to match RNAseq
#'
#' @param target A matrix
#' @param obj A matrix
#' @return A list including the ref and mix normalized and a ggplot
quant_norm_to_target = function(obj, target){
	writeLines("Quantile norm to target..")
	
	if(max(ref$value, na.rm = T) < 50 | max(mix$value, na.rm = T) < 50) stop("ARMET: Both objects have to be in real scale, log tranformation detected")
	
	list(
		ref = target,
		mix = obj %>%
			dplyr:::group_by(sample) %>%
			dplyr:::mutate(
				value = exp(preprocessCore:::normalize.quantiles.use.target(
					log(value+1), 
					target=log(target$value+1)
				))-1
			) %>%
			dplyr:::ungroup(),
		plot = plot(1,1)
	)
	
}

quant_norm_to_RNAseq = quant_norm_to_array = quant_norm_to_target



#' Wrapper function for all normalization based on the current needs
#'
#' @param ref A matrix
#' @param mix A matrix
#' @param is_mix_microarray A cool
#' @return Same list of rnaseq_norm_ref_mix, quant_norm_to_RNAseq
wrapper_normalize_mix_ref = function(mix, ref, is_mix_microarray = F){
	if(!is_mix_microarray)
		rnaseq_norm_ref_mix(ref, mix)
	else
		#quant_norm_to_array(mix, ref)
		quant_norm_to_RNAseq(mix, ref)
}

#' Subsample reference to try homogeneity
#'
#' @param ref A matrix
#' @param cell_types A char vector
#' @param dim_mix A real
#' @return A numerical array of the selected samples
subsample_ref = function(ref, cell_types, dim_mix = 10){
	
	dim_ref_cell_type = min(table(cell_types)) #min(max(table(cell_types)), max(dim_mix/2, 10))
	sam = do.call("c", lapply(unique(cell_types), function(u) {
		which_u = which(cell_types == u)
		if(length(which_u)>1)
			sample(
				which(cell_types == u),
				size = min(dim_ref_cell_type, length(which(cell_types == u)))
			)
		else which_u
	}
	))
	
	return ( sam )
}

#' Plots MDSplot from limma
#'
#' @param mix A matrix
#' @param ref A matrix
#' @param cell_types A char vector
mds_plot = function(mix, ref, cell_types){
	
	do_log_mix =                     if(max(mix, na.rm=T) > 50) T else F
	do_log_ref =                     if(max(ref, na.rm=T) > 50) T else F
	
	if(do_log_mix) mix =             log(mix +1)
	if(do_log_ref) ref =             log(ref+1)
	
	plotMDS(
		cbind(ref, mix),
		labels =                        c(as.character(cell_types), rep("o", dim(mix)[2])),
		col =                           c(rep("black", size=length(cell_types)), rep("red", size=dim(mix)[2]))
	)
	
}

#' Create object for full bayesian model, with the reference samples 
#'
#' @param ref A matrix
#' @param cell_types A char rvector
#' @param markers A char vector
#' @return A list including a two objects
prepare_input = function(ref, node){

	# Sanity check
	error_if_log_transformed(ref)
	
	# Check if sd == 0 for any marker
	ref = check_if_sd_zero_and_correct(ref, node)
	
	# Check if sd == 0 for any marker
	e.obj = ref %>%
		dplyr:::left_join(
			ref %>% 
				distinct(gene) %>% 
				mutate(gene_num = as.numeric(gene)) %>%
				mutate(gene_num = order(gene_num)),
			by = "gene"
		) %>%
		dplyr:::left_join(
			ref %>% 
				distinct(ct) %>% 
				mutate(ct_num = as.numeric(ct)) %>%
				mutate(ct_num = order(ct_num)),
			by = "ct"
		)

	e_mu.obj = e.obj %>%
		dplyr:::distinct(gene_num, ct_num)
	
	# Match e obj with mean e obj
	e.obj = e.obj %>% 
		mutate(
			map_to_mu = 
				match(
					interaction( gene_num, ct_num),
					interaction( e_mu.obj$gene_num, e_mu.obj$ct_num)
				)
		)
		
	list(e=e.obj, e_mu=e_mu.obj)
	
}

#' Averages the signatures of the same cell type
#'
#' @param df A matrix
#' @param do_log A bool
#' @param verbose A bool
#' @return A averaged matrix
get_mean_signature = function(df, do_log=F, verbose= T){
	
	if(nrow(df)==0) stop(sprintf("ARMET: in get_mean_signature the data set has not records for cell types %s. Check bug upstream", paste(colnames(df), collapse=" ")))
	if(length(dim(df))==0) stop(sprintf("ARMET: There are not markers for %s", colnames(df)))
	
	if(max(df, na.rm=T)>50) do_log = T
	if(do_log) df = log(df+1)
	temp =                     data.frame(ct = colnames(df), t(df), row.names=NULL)
	
	
	df.mean =                  do.call("cbind", by(temp, temp$ct, function(df) {
		data.frame(matrixStats:::colMedians(as.matrix(df[,-1, drop=F]), na.rm=T))
	}))
	rownames(df.mean) =        rownames(df)
	
	if(do_log) df.mean = exp(df.mean)-1
	
	# Fix if any cell type is all NA for a gene
	which_with_NA = apply(df.mean, 1, function(mr) any(is.na(mr)))
	if(any(which_with_NA)){
		if(verbose) {
			warning(sprintf("ARMET: in calculating gene means, there are NAs for genes %s\n", paste(rownames(df.mean)[which_with_NA])))
			#	print(df.mean[which_with_NA,])
		} 
		df.mean[is.na(df.mean)] = 0
	} 
	
	return(df.mean)
}

#' Check if prediction of marker value in full ayesian model is too different from the provided reference marker values distibutions
#'
#' @param fit A stan fit object
#' @param map A matrix
check_markers = function(fit, map){
	
	mark =                            summary(fit, pars="x_marker")$summary[,1:2]
	bg =                              summary(fit, pars="b_mu_multiplier")$summary[,c(1,3)]
	e_mu_mu =                      summary(fit, pars="e_mu_mu")$summary[,1]
	e_mu_sigma =                      summary(fit, pars="e_mu_sigma")$summary[,1]
	
	mark.p = apply(mark, 1, function(m){
		# https://stats.stackexchange.com/questions/110225/two-sample-t-test-for-equal-means-with-unequal-variances-for-large-samples
		pnorm(e_mu_mu, m[1], sqrt( e_mu_sigma^2 + m[2]^2 ), lower.tail=F)
	})
	
	bg.p = apply(bg, 1, function(b){
		pnorm(0.8, b[1], b[2], lower.tail=F)
		# library(actuar)
		# y <- rbeta(1000,2,2)
		# loglik <- function(mu, x) { sum(-dbeta(x,mu[1],mu[2],log = TRUE)) }
		# optim(par = c(1,1), fn=loglik,x=y,method = "L-BFGS-B",lower=c(0,0))
	})
	
	# Select bad markers for each cell type
	bad_markers = sapply(colnames(map), function(cn) {
		mc = map[,cn]
		names(mc)[mc==0 & (mark.p<0.05 | bg.p>0.05)]
	})
	
	# Print warning message
	dummy = sapply(1:length(bad_markers), function(i){
		if(length(bad_markers[[i]])>0)
			writeLines(sprintf(
				"ALERT! For cell type %s the marker genes %s is not a stable gene in this data set. Please reconsider.",
				names(bad_markers)[i], paste(unlist(bad_markers[i]), collapse = ", ")
			))
	})
	
}

#' Parse fit object in a 2D matrix
#'
#' @param fit A stan fit object
#' @param param A char
#' @param column A real
parse_result_2D <- function(fit, param, column=1){
	f =                                 summary(fit, param)$summary[,column]
	parse_summary_vector_in_2D          (f)
}

#' Parse a vector from stan fit object in a 2D matrix
#'
#' @param f A array
#' @return A data frame
parse_summary_vector_in_2D = function(f){
	doc.topic =                         do.call("rbind",lapply(names(f),function(x) strsplit(x, "\\[|\\]|,")[[1]]))
	f =
		data.frame(
			topic =                         as.numeric(doc.topic[, 2]),
			gene =                          as.numeric(doc.topic[, 3]),
			mean =                          f
		)
	
	f =                                 reshape:::cast(f, gene ~ topic, value="mean")
	rownames(f) =                       f$gene
	f =                                 data.frame(f[, -1])
	
	return(f)
}

get_stats_on_ref = function(ref, tree){
	`%>%` <- magrittr::`%>%`
	rbind(
		foreach:::foreach(ct = get_leave_names(tree), .combine = rbind) %do% {
			tibble:::tibble(
				ct, 
				count = length(which(get_leave_label(node_from_name(tree, ct), label = "markers") %in%	ref$gene	)),
				val = "gene"
			)
		},
		tibble:::as_tibble(table(ref$ct)) %>%  
			dplyr:::rename(ct=Var1, count = n) %>%   
			dplyr:::mutate(val = "sample")
	)
	
}

plot_model_results = function(){
	
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
	
}