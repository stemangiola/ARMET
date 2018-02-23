
build_data_directory = function(){
	
	# Load reference
	# ref_seq = as.matrix(read.csv("~/PhD/deconvolution/ARMET_dev/ARMET_TME_signature_df_RNAseq.csv", header=T, row.names=1))
	# colnames(ref_seq) = as.vector(sapply(colnames(ref_seq), function(cn) strsplit(cn, ".", fixed=T)[[1]][1]))
	# dic = data.frame(sprintf("s%s", 1:ncol(ref_seq)), colnames(ref_seq))
	# colnames(ref_seq) = dic[,1]
	# ref_seq = tibble::as_tibble(reshape::melt(ref_seq)) %>% dplyr::rename(gene=X1, sample=X2, value = value) %>% dplyr::mutate(ct=dic[match(sample, dic[,1]),2])
	# save(ref_seq, file="data/ref_seq.rda")
	# 
	# ref_array = as.matrix(read.csv("~/PhD/deconvolution/ARMET_dev/ARMET_TME_signature_df_array.csv", header=T, row.names=1))
	# colnames(ref_array) = as.vector(sapply(colnames(ref_array), function(cn) strsplit(cn, ".", fixed=T)[[1]][1]))
	# dic = data.frame(sprintf("s%s", 1:ncol(ref_array)), colnames(ref_array))
	# colnames(ref_array) = dic[,1]
	# ref_array = tibble::as_tibble(reshape::melt(ref_array)) %>% dplyr::rename(gene=X1, sample=X2, value = value) %>% dplyr::mutate(ct=dic[match(sample, dic[,1]),2])
	# save(ref_array, file="data/ref_array.rda")
	
	
	# Create trees
	# library(jsonlite)
	# tree = read_json("/wehisan/home/allstaff/m/mangiola.s/PhD/deconvolution/ARMET_dev/ARMET_TME_tree_RNAseq.json")
	# save(tree, file="data/tree_json.rda")
	
}

average_duplicated_genes_tibble_spreaded = function(tbl){
	
	dup_genes =
		tbl %>%
		dplyr::group_by(gene) %>%
		dplyr::summarise(tot = n()) %>%
		dplyr::filter(tot > 1) %>%
		dplyr::pull(gene) %>%
		as.character()
	
	tbl %>%
		dplyr::mutate(gene = as.character(gene)) %>%
		dplyr::filter(!gene %in% dup_genes) %>%
		dplyr::bind_rows(
			tbl %>%
				dplyr::mutate(gene = as.character(gene)) %>%
				dplyr::filter(gene %in% dup_genes) %>%
				dplyr::group_by(gene) %>%
				dplyr::summarise_if(is.numeric, median) %>%
				dplyr::ungroup()
		)
	
}

#' Check the input for anomalies
check_input = function(mix, is_mix_microarray, my_design, cov_to_test, prior_sd, custom_ref, tree){
	
	if(!tibble::is_tibble(mix)) stop("ARMET: The mixture must be a tibble")
	
	# Check if custom ref is tibble
	if(!is.null(custom_ref)){
		if(!tibble::is_tibble(custom_ref)) 
			stop("ARMET: The reference must be a tibble")
		if(!all(colnames(ref_seq)%in%colnames(custom_ref))) 
			stop(
				sprinf(
					"ARMET: The columns of the provided reference must be %s", 
					paste(colnames(ref_seq), collapse=" ")
				)
			)
	} 
	
	# Check if microarray data
	if(is_mix_microarray) writeLines("ARMET: The input matrix is from microarray.")
	
	# Check if NA in mix
	if(
		mix %>%	
		dplyr::select_if(function(.) any(is.na(.))) %>% 
		ncol() > 0
	) stop("ARMET: NAs found in the query matrix")
	
	# Check if design is tibble
	if(!is.null(cov_to_test) && !tibble::is_tibble(my_design)) stop("ARMET: The design matrix must be a tibble")
	
	# This is how many conditions are in the study (e.g., treatment-vs-non-treatment)
	if(!is.null(my_design)) {
		writeLines("ARMET: The design matrix is :")
		print(head(my_design))
	}
	
	# Set up the covariate to test if any
	if(is.null(my_design) & !is.null(cov_to_test)) stop("ARMET: you have specified a covariate to test but no design matrix")
	
	# Check prior sd
	if(prior_sd<=0) stop("ARMET: prior_sd must be a positive number.")
	
	
	
	# Check custom reference
	if(
		!is.null(custom_ref) &&
		!all(get_leave_names(tree, last_level = 1) %in% (custom_ref %>% dplyr::distinct(ct)  %>% dplyr::pull(ct) %>% as.character()))
	) {
		writeLines("ARMET: Some cell types within the tree are absent in your custom reference provided") 
		cts = get_leave_names(tree, last_level = 1) 
		print( cts[!cts %in% (custom_ref %>% dplyr::distinct(ct) %>% dplyr::pull(ct) %>% as.character())	]	)
		
		stop()
	}
	
	# Check if duplicated genes in mix
	if(
		mix %>%
		dplyr::group_by(gene) %>%
		dplyr::summarise(tot = n()) %>%
		dplyr::filter(tot > 1) %>%
		nrow() > 1
	) {
		writeLines("ARMET: There are duplicated genes n the provided mix. Genes will be summarised with medians.")
		
		pf = parent.frame()
		
		pf$mix = average_duplicated_genes_tibble_spreaded(mix)
	}
	
	
}

#' Get parameters from ini file if it exists
get_ini = function(){
	if(file.exists("ARMET.ini")) {
		pars =  ini::read.ini("ARMET.ini")$ARMET
		writeLines("ARMET: Importing parameters from .ini file..")
		for(n in names(pars)){
			if(pars[n]%in%c("T", "F", "TRUE", "FALSE"))
				pars[n] = as.logical(pars[n])
			else 	pars[n] = strsplit(as.character(pars[n]), ",")
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
	
	cl = parallel::makeCluster(20)
	parallel::clusterExport(cl, "temp", envir = environment())
	parallel::clusterEvalQ(cl, library("matrixStats"))
	
	ref.symbol = do.call("rbind", parallel::parLapply(cl, us, function(tt) {
		matrixStats::colMedians(as.matrix(temp[temp$symbol==tt,-1]))
	}))
	parallel::stopCluster(cl)
	
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
	reference.4plot = data.frame(reshape::melt(reference.4plot), source="reference")
	obj.4plot = obj
	colnames(obj.4plot) = paste0("r", 1:dim(obj.4plot)[2])
	obj.4plot = data.frame(reshape::melt(obj.4plot), source="obj")
	
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
	
	plot_chains = rstan::traceplot(fit, inc_warmup = FALSE, pars="alpha")
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
		dplyr::group_by(gene, ct) %>%
		dplyr::summarise(
			log_sigma = stats::mad(log(value+1)), 
			log_avg = mean(log(value+1))
		) %>%
		dplyr::ungroup() %>%
		dplyr::mutate(to_recalculate = log_sigma==0) 
	
	find_closest_sd_uo_the_tree = function(my_ct, my_gene, tb){
		#	print(my_ct)
		#	print(my_gene)
		
		# Iterate upwards thoward root to get the sd from the closest group
		hierarchy = get_hierarchy(node, my_ct)
		ancestors = rev(hierarchy[1:(length(hierarchy)-1)])
		
		for(h in ancestors){
			
			# Get first descendant of ancestor including the ancestor 
			ct_to_consider = unlist(c(h, get_leave_label(node_from_name(node, h), recursive = F)))
			
			mean_log_sigma = tb %>% 
				dplyr::filter(
					ct %in% ct_to_consider & 
						gene == my_gene & 
						log_sigma > 0) %>%
				
						{
							
							# Exception if I have only one sample per cell type
							if((.) %>% nrow() == 0) 
								(.) %>% dplyr::summarise(mean_log_sigma = 0)
							else 
								(.) %>% dplyr::summarise(mean_log_sigma = mean(log_sigma))
							
						}
			
			if(mean_log_sigma>0) break
			
		}
		
		# If I have only one sample per cell type output 1 as SD
		ifelse(mean_log_sigma == 0, 1, as.numeric(mean_log_sigma))
		
	}
	
	# Add closest sd
	df.summary = df.summary %>%
		dplyr::rowwise() %>%
		dplyr::mutate(
			log_sigma = 
				ifelse(
					log_sigma==0,
					find_closest_sd_uo_the_tree(ct, gene, df.summary),
					log_sigma
				)
		) %>%
		dplyr::ungroup()
	
	if(nrow(df.summary %>% dplyr::filter(to_recalculate))>0){
		writeLines("ARMET: One of your markers has 0 variance, this is not compatible with MCMC inference.")
		writeLines("ARMET: The following genes will acquire a background variance.")
		print( df.summary %>% dplyr::filter(to_recalculate) )
	}
	
	# Sample for sd == 0
	df %>%
		dplyr::left_join(df.summary, by=c("gene", "ct")) %>%
		dplyr::rowwise() %>%
		dplyr::mutate(value = ifelse(to_recalculate, exp(rnorm(1, log_avg, log_sigma)), value)) %>%
		dplyr::ungroup() %>%
		dplyr::select(-log_sigma, -log_avg, -to_recalculate)
	
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
	{ if("ct"%in%names(.)) dplyr::select(-ct)	else .}	%>%
		tidyr::spread(sample, value) %>% 
		tidyr::drop_na() %>%
		dplyr::do(
			(.) %>% 
				dplyr::filter(
					rowSums(
						edgeR::cpm(
							(.) %>% 
								dplyr::select(-gene)
						) > cpm_theshold
					) >= 
						ceiling(ncol( 
							(.) %>% dplyr::select(-gene)
						) * prop)
				)
		)
	
	list(
		nf = tibble::tibble(
			sample = factor(colnames(df.filt %>% dplyr::select(-gene))),
			nf = edgeR::calcNormFactors(df.filt %>% dplyr::select(-gene), refColumn=reference)
		),
		df = df.filt %>% 
			tidyr::gather(sample, value, 2:ncol(.)) %>% 
			dplyr::mutate(sample = factor(sample)) %>%
			{ 
				if("ct"%in%names(.)) 
					dplyr::left_join(
						df %>% dplyr::distinct(sample, ct), 
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
	#if(verbose) writeLines("Normalizing RNA-seq data with TMM")
	
	nf.obj = rnaseq_norm.calcNormFactor(df, reference, cpm_theshold, prop)
	nf = nf.obj$nf
	df.filt = nf.obj$df
	
	df %>% 
		dplyr::group_by(sample) %>% 
		dplyr::mutate(tot = sum(value, na.rm = T)) %>%
		dplyr::ungroup() %>%
		dplyr::left_join(nf, by="sample") %>%
		dplyr::left_join(
			df.filt %>%
				dplyr::group_by(sample) %>%
				dplyr::summarise(tot_filt = sum(value, na.rm = T)),
			"sample"
		) %>%
		dplyr::mutate(value = value / (tot_filt * nf) * max(tot)) %>%
		dplyr::select(-tot, -nf, -tot_filt)
	
}

#' Normalize ref to match mix using TMM
#'
#' @param ref A matrix
#' @param mix A matrix
#' @return A list including the ref and mix normalized and a ggplot
rnaseq_norm_ref_mix = function(obj, target){
	
	error_if_log_transformed(target)
	error_if_log_transformed(obj)
	
	# Normalize the obj
	obj = rnaseq_norm(obj)
	
	# Calclate normalization factors for obj and target => two numbers
	nf.obj = rnaseq_norm.calcNormFactor(
		rbind(
			obj %>% 
				dplyr::group_by(gene) %>% 
				dplyr::summarise(value = median(value, na.rm=T)) %>% 
				dplyr::mutate(sample="obj"),
			target %>% 
				dplyr::group_by(gene) %>% 
				dplyr::summarise(value = median(value, na.rm=T)) %>% 
				dplyr::mutate(sample="target")
		), 
		1)
	
	nf = nf.obj$nf
	my_df = nf.obj$df
	
	tot_reference = (my_df %>%
									 	dplyr::filter(sample == "obj") %>% 
									 	dplyr::summarise(tot = sum(value, na.rm=T)))$tot 
	
	tot_other = (my_df %>%
							 	dplyr::filter(sample == "target") %>% 
							 	dplyr::summarise(tot = sum(value, na.rm=T))
	)$tot 
	
	obj = obj %>%
		dplyr::mutate(tot_ref = tot_reference) %>%
		dplyr::mutate(nf = 
										(nf %>%
										 	dplyr::filter(sample=="obj"))$nf
		) %>%
		dplyr::mutate(tot = tot_ref) %>%
		dplyr::mutate(value = value / (tot * nf) * tot_ref) %>%
		dplyr::select(-tot, -tot_ref, -nf)
	
	target = target %>%
		dplyr::mutate(tot_ref = tot_reference) %>%
		dplyr::mutate(nf = 
										(nf %>%
										 	dplyr::filter(sample=="target"))$nf
		) %>%
		dplyr::mutate(tot = tot_other) %>%
		dplyr::mutate(value = value / (tot * nf) * tot_ref) %>%
		dplyr::select(-tot, -tot_ref, -nf)
	
	p = plot_densities( 
		rbind( 
			target %>% 
				dplyr::mutate(color=1) %>%
				dplyr::select(-ct),
			obj %>% 
				dplyr::mutate(color=1)
		)
	)
	
	return(list(target = target, obj = obj,  plot = p))
}

#' Normalize array to match RNAseq
#'
#' @param target A matrix
#' @param obj A matrix
#' @return A list including the ref and mix normalized and a ggplot
quant_norm_to_target = function(obj, target){
	#writeLines("ARMET: Quantile normalization")
	
	error_if_log_transformed(obj)
	error_if_log_transformed(target)
	
	# Needed for: Transform tibble in matrix and normalize because is faster
	obj = obj %>%
		tidyr::spread(sample, value)
	
	obj = dplyr::bind_cols(
		obj %>% dplyr::select(gene),
		
		data.frame(
			exp(preprocessCore::normalize.quantiles.use.target(
				as.matrix(log(obj[,-1]+1)), 
				target=log(target$value+1)
			)) - 1
		) %>% 
			tibble::as_tibble()
	) %>%
		magrittr::set_colnames(colnames(obj)) %>%
		tidyr::gather(sample, value, -gene) %>%
		dplyr::mutate_if(is.character, as.factor)
	
	list(
		target = target,
		obj = obj,
		plot = plot(1,1)
	)
	
}

#' Wrapper function for all normalization based on the current needs
#'
#' @param ref A matrix
#' @param mix A matrix
#' @param is_mix_microarray A cool
#' @return Same list of rnaseq_norm_ref_mix, quant_norm_to_target
wrapper_normalize_mix_ref = function(mix, ref, is_mix_microarray){
	if(!is_mix_microarray) {
		norm = rnaseq_norm_ref_mix(mix, ref)
		list(mix=norm$obj, ref=norm$target)
	}
	else {
		norm = quant_norm_to_target(mix, ref)
		list(mix=norm$obj, ref=norm$target)
	}
	
	
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
		dplyr::left_join(
			ref %>% 
				dplyr::distinct(gene) %>% 
				dplyr::mutate(gene_num = as.numeric(gene)) %>%
				dplyr::mutate(gene_num = order(gene_num)),
			by = "gene"
		) %>%
		dplyr::left_join(
			ref %>% 
				dplyr::distinct(ct) %>% 
				dplyr::mutate(ct_num = as.numeric(ct)) %>%
				dplyr::mutate(ct_num = order(ct_num)),
			by = "ct"
		)
	
	e_mu.obj = e.obj %>%
		dplyr::distinct(gene_num, ct_num)
	
	# Match e obj with mean e obj
	e.obj = e.obj %>% 
		dplyr::mutate(
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
		data.frame(matrixStats::colMedians(as.matrix(df[,-1, drop=F]), na.rm=T))
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
	
	f =                                 reshape::cast(f, gene ~ topic, value="mean")
	rownames(f) =                       f$gene
	f =                                 data.frame(f[, -1])
	
	return(f)
}

get_stats_on_ref = function(ref, tree){
	
	rbind(
		foreach::foreach(ct = get_leave_names(tree), .combine = rbind) %do% {
			tibble::tibble(
				ct, 
				count = length(which(get_leave_label(node_from_name(tree, ct), label = "markers") %in%	ref$gene	)),
				val = "gene"
			)
		},
		tibble::as_tibble(table(ref$ct)) %>%  
			dplyr::rename(ct=Var1, count = n) %>%   
			dplyr::mutate(val = "sample")
	)
	
}

plot_model_results = function(){
	
	proportions_4_plot = as.data.frame(proportions)
	proportions_4_plot$sample = rownames(proportions_4_plot)
	proportions_4_plot$cov = if(length(cov_to_test)>0) my_design[,"cov"] else 1
	proportions_4_plot = reshape::melt(proportions_4_plot, id.vars=c("sample", "cov"))
	proportions_4_plot$perc_1 = reshape::melt(proportions_2.5)$value
	proportions_4_plot$perc_2 = reshape::melt(proportions_97.5)$value
	proportions_4_plot$minus_sd = proportions_4_plot$value - reshape::melt(proportions_sd)$value
	proportions_4_plot$plus_sd = proportions_4_plot$value + reshape::melt(proportions_sd)$value
	p = ggplot2::ggplot( proportions_4_plot, ggplot2::aes(x=jitter(cov), y=value,fill=factor(variable)))+ 
		ggplot2::geom_boxplot(coef = 6) + 
		ggplot2::geom_point(size=0.1) +
		ggplot2::geom_linerange(ggplot2::aes(ymin = perc_1, ymax = perc_2),  alpha=0.05) +
		ggplot2::geom_errorbar(ggplot2::aes(ymin = minus_sd, ymax = plus_sd),width=.2,  alpha=0.2) +
		ggplot2::facet_grid(~variable) +
		ggplot2::theme(
			axis.text.x= ggplot2::element_text(angle=90, vjust=0.4,hjust=1), 
			panel.background = ggplot2::element_blank()
		)
	
	#if(save_report) ggplot2::ggsave(sprintf("%s/%s_proportions.pdf", output_dir, ct), useDingbats=FALSE,plot=p)
	
	if(length(cov_to_test)>0){
		#plot generate quantities
		gq <- as.matrix(fit, pars = c("beta_gen"))[1,]
		gq= parse_summary_vector_in_2D(gq)
		colnames(gq) = colnames(mix)
		rownames(gq) = colnames(model.in$x)
		gq = as.data.frame(t(gq))
		gq$cov = my_design[,"cov"]
		gq$sample = rownames(gq)
		gq = reshape::melt(gq, id.vars=c("sample", "cov"))
		p = ggplot2::ggplot( gq, ggplot2::aes(x=factor(cov), y=value,fill=factor(variable)))+
			ggplot2::geom_boxplot()+ 
			ggplot2::geom_jitter(position=position_dodge(width=0.75), ggplot2::aes(group=variable)) +
			ggplot2::facet_grid(~variable) + 
			ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, vjust=0.4,hjust=1))
		
		#if(save_report) ggplot2::ggsave(sprintf("%s/%s_proportions_posterior.pdf", output_dir, ct), useDingbats=FALSE,plot=p)
	}
	# Calculate predicted values
	#pred.df=as.matrix(ref[markers,])%*%t(as.matrix(res.df))
	
	
	# about means
	df_4_plot = data.frame(
		merge(
			reshape::melt(mix+1), 
			reshape::melt(t(proportions %*% t(ref.mean)+1)), 
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
	
	p = ggplot2::ggplot(df_4_plot, ggplot2::aes(predicted, value, label = X1, size=bg_prop,color = factor(ct))) + 
		ggplot2::geom_point(alpha=0.5) + 
		ggplot2::scale_size(range = c(0, 5)) + 
		ggplot2::expand_limits(x = 0.1, y = 0.1) +
		ggplot2::geom_abline(intercept=0, slope=1) + 
		ggplot2::scale_x_log10() + 
		ggplot2::scale_y_log10() +
		ggrepel::geom_text_repel(
			ggplot2::aes(color = ct),
			size = 2,
			segment.alpha= 0.2) +
		ggplot2::facet_wrap( ~ X2, nrow=5) +
		ggplot2::coord_fixed() +
		ggplot2::theme(
			panel.background = ggplot2::element_blank(), 
			legend.text=ggplot2::element_text(size=30)
		)
	
	#if(save_report) ggplot2::ggsave(sprintf("%s/%s_predicted_values.pdf", output_dir, ct), useDingbats=FALSE, width = 80, height = 80, units = "cm", plot=p)
	
	
	if(!is.null(observed_prop) & 0) {
		writeLines("ARMET: start plotting expected gene vallues")
		link = tibble::as_tibble(do.call("rbind", lapply(unique(colnames(observed_prop)), function(op) {
			ct = get_hierarchy(my_tree, ct=op)
			c(op, ct[ct%in%colnames(ref.mean)])
		})))
		if(ncol(link)>=2){
			colnames(link) = c("ct", "ct_new")
			observed_prop = tibble::as_tibble(reshape::melt(observed_prop))
			colnames(observed_prop) = c("sample", "ct", "value")
			
			observed_prop = 
				dplyr::left_join(observed_prop, link, by="ct") %>% 
				dplyr::select(-ct) %>% 
				dplyr::group_by(sample, ct_new) %>% 
				dplyr::summarise(value = sum(value, na.rm=TRUE)) %>% 
				tidyr::spread(ct_new, value)
			
			for(r in colnames(ref.mean)[!colnames(ref.mean)%in%colnames(observed_prop)]) {
				my_df = matrix(rep(0, nrow(observed_prop)))
				colnames(my_df) = r
				observed_prop = observed_prop %>% dplyr::bind_cols(as_tibble(my_df)) 
			}
			
			observed_prop = as.data.frame(observed_prop)
			rownames(observed_prop) = observed_prop[,1]
			observed_prop = observed_prop[,-1]
			observed_prop = as.matrix(observed_prop)
			observed_prop = observed_prop[,colnames(ref.mean)]
			
			print(proportions)
			
			df_4_plot = data.frame(merge(reshape::melt(mix+1), reshape::melt(t(observed_prop %*% t(ref.mean)+1)), by=c("X1", "X2")))
			colnames(df_4_plot)[3:4] = c("value", "predicted")
			df_4_plot$ct = sapply(df_4_plot$X1, function(x1) node_from_gene(my_tree, x1)$name )
			df_4_plot$bg = apply(df_4_plot, 1, function(mr) 	y_hat_background[mr[2], mr[1]] )
			df_4_plot$bg_prop = df_4_plot$bg /df_4_plot$value
			df_4_plot$bg_prop[df_4_plot$bg_prop >10] = 10
			
			p = ggplot2::ggplot( df_4_plot[df_4_plot$X2%in%unique(df_4_plot$X2)[1:20], ], ggplot2::aes(predicted, value, label = X1, size=bg_prop)) + 
				ggplot2::geom_point(alpha=0.5) + 
				ggplot2::scale_size(range = c(0, 5)) + 
				ggplot2::expand_limits(x = 0.1, y = 0.1) +
				ggplot2::geom_abline(intercept=0, slope=1) + 
				ggplot2::scale_x_log10() + 
				ggplot2::scale_y_log10() +
				ggrepel::geom_text_repel(
					aes(color = ct),
					size = 2,
					segment.alpha= 0.2) +
				ggplot2::facet_wrap( ~ X2, nrow=5) +
				ggplot2::coord_fixed() +
				ggplot2::theme(
					panel.background = ggplot2::element_blank(), 
					legend.text=ggplot2::element_text(size=30)
				)
			
			#if(save_report) ggplot2::ggsave(sprintf("%s/%s_predicted_values_observed_prop.pdf", output_dir, ct), useDingbats=FALSE, width = 80, height = 80, units = "cm", plot=p)
		}
		
	}
	
}



get_center_bg = function(coef_ang_posterior){
	
	my_sd = matrixStats::colSds(coef_ang_posterior)
	names(my_sd) = colnames(coef_ang_posterior)
	my_mean = colMeans(coef_ang_posterior)
	if(length(my_mean) <= 2) return(list(mean =mean(my_mean), sd= sum(my_sd)))
	
	my_dist = stats::dist(my_mean, diag = T, upper = T)
	my_dist.m = as.matrix(my_dist)
	my_clust = stats::hclust(my_dist)
	
	w = which(my_dist.m==min(my_dist.m[my_dist.m>0])  )[1]
	my_row = floor(w/ncol(my_dist.m))
	if(w-(my_row*ncol(my_dist.m)) == 0) {
		my_col =  ncol(my_dist.m)
	} else {
		my_col =  w - my_row*ncol(my_dist.m)
		my_row = my_row + 1
	}
	
	my_couple = as.data.frame(my_dist.m[my_row, my_col, drop=F])
	my_couple = c(rownames(my_couple), colnames(my_couple))
	
	# Find value of the supposely background standard deviation
	my_hierarchical_sd = sum(my_sd[my_couple])
	
	# Cue the tree for cluster closer than 95 percentile of that standard deviation
	my_cut = stats::cutree(my_clust, h=my_hierarchical_sd*2)
	
	# Calculate cluster size
	my_table = table(my_cut)
	#print(my_table)
	
	# Check if anything clustered at all
	
	if(any(my_table>1)) {
		
		# Get bg cluster among possibly many
		my_c = unique(my_cut[names(my_cut)%in%my_couple])
		
		# Get elements in that cluster
		my_elem = names(my_cut)[my_cut==my_c]
		
		if(length(my_c)>1) stop("ARMET: the closest samples do not form cluster even though a cluster is present")
		
		my_new_table = my_table
		my_table_check = NULL
		while(any(my_new_table>1) & length(my_new_table)>2 & any(my_table_check != my_new_table)){
			my_table_check = my_new_table
			my_ancestor = mean(my_mean[my_couple])
			names(my_ancestor) = paste(my_couple, collapse="_")
			my_new_mean = c(my_mean[!names(my_mean)%in%my_couple], my_ancestor)
			my_new_dist = stats::dist(my_new_mean, diag = T, upper = T)
			my_new_dist.m = as.matrix(my_new_dist)
			my_new_clust = stats::hclust(my_new_dist)
			my_new_cut = stats::cutree(my_new_clust, h=my_hierarchical_sd*2)
			my_new_c = my_new_cut[names(my_new_cut)%in%names(my_ancestor)]
			my_new_elem = names(my_new_cut)[my_new_cut==my_new_c]
			my_elem = c(my_elem, my_new_elem[!my_new_elem%in%names(my_ancestor)])
			my_new_table = table(my_new_cut)
			my_hierarchical_sd = sum(my_sd[my_elem])
		}
		list(mean = mean(my_mean[my_elem]), sd= my_hierarchical_sd)
	}
	else list(mean = median(my_mean), sd = my_sd[my_couple])
	
}

#' Hipotesis test for the covariate of choice
#' @rdname dirReg_test
#'
#' Prints a report of the hipothesis testing
#'
#' @param fit stan fit object 
#' @param my_design design matrix
#' @param cov_to_test character variable indicating the column name of the design matrix to test
#' @param which_cov_to_test for internal usage
#'
#' @return a vector including
#'     pvalues of the intercepts
#'     pvalues of the angular coefficients
#'     sign of the angular coefficient
#'
#' @examples
#'  dirReg_test(fit, my_design, cov_to_test)
#' @export
dirReg_test = function(fit, my_design, cov_to_test = NULL, which_cov_to_test = 2, names_groups = NULL){

	logit_adj <- function(v, t=0.5) -log(t*(v-1) / ((t-1)*v));
	
	# Decide which covariate to check
	if(!is.null(cov_to_test)) which_cov_to_test = which(colnames( my_design %>% dplyr::select(-sample) ) == cov_to_test )
	
	# Get the posterior distributions
	alpha_posterior = as.matrix(as.data.frame(rstan:::extract( fit, "alpha")))
	coef_ang_posterior = as.matrix(alpha_posterior[,grep(sprintf("alpha.%s", which_cov_to_test), colnames(alpha_posterior), fixed=T)])
	interc_posterior = as.matrix(alpha_posterior[,grep("alpha.1", colnames(alpha_posterior), fixed=T)])
	
	K = ncol(interc_posterior)
	gcb = get_center_bg(coef_ang_posterior)
	m = gcb$mean
	
	# plot = plot_densities(logit_adj(coef_ang_posterior, m)  , do_log = F, color="0",  fill = 1:K) +
	# 	ggtitle("Causal trends - posterior distribution of angular coeff. - Simplex regression")
	# 
	stats = do.call("rbind", lapply(1:K, function(i){
		
		mcap = mean(coef_ang_posterior[,i])
		scap = sd(coef_ang_posterior[,i])
		ecap = scap/sqrt(nrow(alpha_posterior))
		mip = mean(interc_posterior[,i])
		sip = sd(interc_posterior[,i])
		eip = sip/sqrt(nrow(alpha_posterior))
		pcap = min( 2*(1-pnorm(m, mcap, scap)), 2*(1-pnorm(m, mcap, scap, lower.tail = F)))
		pip = min ( 2*(1-pnorm(1/K, mip, sip)), 2*(1-pnorm(1/K, mip, sip, lower.tail = F)))
		
		pcap = min( 2*(1-pnorm(mcap, gcb$mean, gcb$sd)), 2*(1-pnorm(mcap, gcb$mean, gcb$sd, lower.tail = F)))
		
		mip_human_readable = mip
		mcap_human_readable = logit_adj(mcap, m)
		pcap_human_readable = if(pcap<0.01 & pcap>0) formatC(pcap, format = "e", digits = 2) else round(pcap, 2)
		pip_human_readable = if(pip<0.01& pip>0) formatC(pip, format = "e", digits = 2) else round(pip, 2)
		
		tibble::tibble(
			m = m,
			mcap = mcap,
			scap = scap,
			ecap = ecap,
			mip = mip,
			sip = sip,
			eip = eip,
			pcap = pcap,
			pip = pip,
			mip_human_readable = mip_human_readable,
			mcap_human_readable = mcap_human_readable,
			pcap_human_readable = pcap_human_readable,
			pip_human_readable = pip_human_readable,
			symbol_pcap = if(pcap<0.001) "***" else if(pcap<0.01) "**" else if(pcap<0.05) "*" else if(pcap<0.1) "." else "",
			symbol_pip = if(pip<0.001) "***" else if(pip<0.01) "**" else if(pip<0.05) "*" else if(pip<0.1) "." else "",
			direction = if(mcap_human_readable>0) "+" else if(mcap_human_readable<0) "-" else "0"
		)
		
	}))
	stats$ct = names_groups
	
	stats = stats %>% dplyr::mutate_if(is.character, as.factor)
	
	list(stats=stats, plot=plot)
}

parse_extract_2D <- function(fit, param, fun)	apply(rstan::extract(fit, param)[[1]], c(2,3), fun )

#' @rdname beta_reg_hierarchical
#' @export
beta_reg_hierarchical = function(fit, my_design){
	browser()
	
	writeLines("ARMET: Starting inference beta reg..")
	
	beta_posterior_mu = parse_extract_2D(fit, "beta", mean)
	beta_posterior_sd = parse_extract_2D(fit, "beta", sd)
	beta_posterior = aperm(rstan::extract(fit, "beta")[[1]], c(2,3,1))
	beta_posterior = beta_posterior[,,1:100]
	
	rstan:::sampling( stanmodels$beta_reg_hierarchical, 
										#stan(file="~/PhD/simplexRegression/src/stan_files/beta_reg_hierarchical.stan",
										data=list(
											I = dim(beta_posterior)[3],
											K=ncol(beta_posterior_mu),
											N=nrow(beta_posterior_mu),
											X=my_design,
											R=ncol(my_design),
											beta_posterior = beta_posterior,
											beta_mu=beta_posterior_mu,
											beta_sd = beta_posterior_sd
										),
										cores=4,
										iter = 1000,
										refresh = 0
	)
}

#' Hipotesis test for the covariate of choice
#' @rdname betaReg_test
#'
#' Prints a report of the hipothesis testing
#'
#' @param fit stan fit object 
#' @param my_design design matrix
#' @param cov_to_test character variable indicating the column name of the design matrix to test
#' @param which_cov_to_test for internal usage
#'
#' @return a vector including
#'     pvalues of the intercepts
#'     pvalues of the angular coefficients
#'     sign of the angular coefficient
#'
#' @examples
#'  betaReg_test(fit, my_design, cov_to_test)
#' @export
betaReg_test = function(fit, my_design, cov_to_test = NULL, which_cov_to_test = 2, names_groups = NULL){
	browser()
	# Decide which covariate to check
	if(!is.null(cov_to_test)) which_cov_to_test = which(colnames( my_design %>% dplyr::select(-sample) ) == cov_to_test )
	
	# Get the posterior distributions
	alpha_posterior = as.data.frame(rstan:::extract( fit, "alpha2"))
	coef_ang_posterior = as.matrix(alpha_posterior[,grep(sprintf("alpha2.%s", which_cov_to_test), colnames(alpha_posterior), fixed=T)])
	interc_posterior = as.matrix(alpha_posterior[,grep("alpha2.1", colnames(alpha_posterior), fixed=T)])
	
	K = ncol(interc_posterior)
	m = 0
	
	plot = plot_densities(coef_ang_posterior , do_log = F, color="0",  fill = 1:K) +
		ggtitle("Effectual trends - posterior distribution of angular coeff. - Beta regression")
	
	stats = 	do.call("rbind", lapply(1:K, function(i){
		
		mcap = mean(coef_ang_posterior[,i])
		scap = sd(coef_ang_posterior[,i])
		ecap = scap/sqrt(nrow(alpha_posterior))
		mip = mean(interc_posterior[,i])
		sip = sd(interc_posterior[,i])
		eip = sip/sqrt(nrow(alpha_posterior))
		pcap = min( 2*(1-pnorm(m, mcap, scap)), 2*(1-pnorm(m, mcap, scap, lower.tail = F)))
		pip = min ( 2*(1-pnorm(1/K, mip, sip)), 2*(1-pnorm(1/K, mip, sip, lower.tail = F)))
		mip_human_readable = inv_logit(mip)
		mcap_human_readable = mcap
		pcap_human_readable = if(pcap<0.01 & pcap>0) formatC(pcap, format = "e", digits = 2) else round(pcap, 2)
		pip_human_readable = if(pip<0.01& pip>0) formatC(pip, format = "e", digits = 2) else round(pip, 2)
		
		tibble(
			m = m,
			mcap = mcap,
			scap = scap,
			ecap = ecap,
			mip = mip,
			sip = sip,
			eip = eip,
			pcap = pcap,
			pip = pip,
			mip_human_readable = mip_human_readable,
			mcap_human_readable = mcap_human_readable,
			pcap_human_readable = pcap_human_readable,
			pip_human_readable = pip_human_readable,
			symbol_pcap = if(pcap<0.001) "***" else if(pcap<0.01) "**" else if(pcap<0.05) "*" else if(pcap<0.1) "." else "",
			symbol_pip = if(pip<0.001) "***" else if(pip<0.01) "**" else if(pip<0.05) "*" else if(pip<0.1) "." else "",
			direction = if(mcap_human_readable>0) "+" else if(mcap_human_readable<0) "-" else "0"
		)
		
	}))
	stats$ct = names_groups
	
	stats = stats %>% dplyr::mutate_if(is.character, as.factor)
	
	list(stats=stats, plot=plot)
	
}