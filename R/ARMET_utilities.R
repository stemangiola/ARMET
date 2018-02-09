
#' Check the input for anomalies
check_input = function(mix, is_mix_microarray, my_design, cov_to_test, prior_sd){
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
	
	if(prior_sd<=0) stop("ARMET: prior_sd must be a positive number.")
	
	
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
	
	df = as.matrix(df)
	df.4plot = df
	if(any(table(colnames(df.4plot))>1) | is.null(colnames(df.4plot))){
		colnames(df.4plot) = paste0("#", 1:dim(df.4plot)[2])
		warning("ARMET: colnames duplicated, replacing with dummy names.")
	}
	color_df = tibble:::tibble(X2 = factor(colnames(df.4plot)), color = color, fill = factor(fill))
	df.4plot.melt = tibble:::as_tibble(reshape:::melt(df.4plot))
	df.4plot.melt = dplyr:::left_join(df.4plot.melt, color_df, by="X2")
	if(do_log) df.4plot.melt$value = df.4plot.melt$value + 0.1
	
	p = ggplot2:::ggplot(df.4plot.melt, ggplot2:::aes(value, group=X2, color = color, fill = fill)) +
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
	
	p = ggplot2:::ggplot(df_melted, aes(value+0.1, group=X2, color=source)) +
		geom_line(stat="density", alpha=0.15) +
		scale_x_log10() +
		expand_limits(x=0.1) +
		theme_bw() +
		theme(
			panel.border = element_blank(),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			axis.line = element_line(colour = "black")
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
	reference.4plot = data.frame(melt(reference.4plot), source="reference")
	obj.4plot = obj
	colnames(obj.4plot) = paste0("r", 1:dim(obj.4plot)[2])
	obj.4plot = data.frame(melt(obj.4plot), source="obj")
	
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

#' Checks if the standard deviation of a gene is 0
#'
#' @param df A matrix
check_if_sd_zero = function(df){
	
	if(is.null(df)) stop("ARMET: check_if_sd_zero: The input is null.. please debug")
	if(max(df, na.rm=T)<50) stop("ARMET: check_if_sd_zero: The input was log transformed")
	
	# check if pathologic data with no sd
	df = as.matrix(df)
	df.melt =                        reshape:::melt(df)
	colnames(df.melt) =              c("gene", "cell_type", "value")
	#df.melt$gene.numeric =           match(df.melt$gene, rownames(df))
	#df.melt$cell_type.numeric =      match(df.melt$cell_type, levels(cell_types))
	df.melt.aggr = stats:::aggregate(log(df.melt$value+1), by=list(df.melt$gene, df.melt$cell_type), FUN=function(x) c(mean = mean(x, na.rm=T), sd = sd(x, na.rm=T) ))
	
	# Fix all NA problem
	if(any(apply(is.na(df.melt.aggr$x[,1:2]), 1, function(mr) mr[1] | mr[2] )))
		df.melt.aggr[apply(is.na(df.melt.aggr$x[,1:2]), 1, function(mr) mr[1] | mr[2] ),]$x[,1:2] = 0
	
	if(any(df.melt.aggr[,3][,2]==0)){
		
		writeLines("ARMET: One of your markers has 0 variance, this is not compatible with MCMC inference.")
		writeLines("The following genes will acquire a background variance.")
		df.melt.aggr.2correct = df.melt.aggr[df.melt.aggr[,3][,2]==0,]
		
		
		
		df.melt.aggr.2correct = t(apply(df.melt.aggr.2correct, 1, function(mr){
			close = subset(df.melt.aggr, Group.1==mr[1])
			close = close[order(abs(close[,3][,1] - as.numeric(mr[3]))),]
			new_sd = mean(close[close[,3][,2]>0,][,3][,2][1:5])
			new_sd = 0.1
			c(mr[1:3],new_sd)
		}))
		
		for(i in 1:dim(df.melt.aggr.2correct)[1]){
			ml = length(df[rownames(df)==df.melt.aggr.2correct[i,1], colnames(df)==df.melt.aggr.2correct[i,2]])
			new_log_mean = as.numeric(df.melt.aggr.2correct[i,3])
			new_log_sd = as.numeric(df.melt.aggr.2correct[i,4])
			#print(c(df.melt.aggr.2correct[i,], new_log_mean = new_log_mean, new_log_sd = new_log_sd))
			df[rownames(df)==df.melt.aggr.2correct[i,1], colnames(df)==df.melt.aggr.2correct[i,2]] = exp(abs(rnorm(ml, new_log_mean , new_log_sd)))
		}
	}
	
	df
	
}

#' Calculate the norm factor with calcNormFactor from limma
#'
#' @param df A matrix
#' @param reference A reference matrix
#' @param cpm_theshold A number
#' @param prop A number
#' @return A list including the filtered data frame and the normalization factors
rnaseq_norm.calcNormFactor = function(df, reference = NULL, cpm_theshold = 0.5, prop = 3/4){
	#print(dim(df))
	
	my_df = na.omit(df)
	#print(dim(my_df))
	if(max(my_df, na.rm = T) < 50) stop("ARMET: Both mixture and signatures have to be in count form, log tranformation detected")
	
	cn = colnames(my_df)
	keep1 <- rowSums(edgeR:::cpm(my_df) > cpm_theshold) >= ceiling(ncol(my_df)*prop)
	writeLines(sprintf("ARMET: %s genes on %s total genes have been filtered out for normalization", length(which(!keep1)), nrow(my_df)))
	#print(dim(my_df))
	#print(length(keep1))
	my_df = my_df[keep1, , drop=FALSE]
	
	list(nf=edgeR:::calcNormFactors(my_df, refColumn=reference), df=my_df)
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
	#print(head(df))
	nf.obj = rnaseq_norm.calcNormFactor(df, reference, cpm_theshold, prop)
	nf = nf.obj$nf
	my_df = nf.obj$df
	
	df.norm =
		do.call(
			"cbind",
			lapply(
				1:dim(df)[2],
				function(i) df[,i] / (sum(my_df[,i], na.rm=T) * nf[i])
			)
		)
	
	if(is.null(reference)) reference_value = max(colSums(df, na.rm=T), na.rm = T)
	else reference_value = sum(reference)
	df.norm = df.norm * reference_value
	colnames(df.norm) = colnames(my_df)
	#p = plot_densities_double(reference,obj)
	
	return(df.norm)
}

#' Normalize ref to match mix using TMM
#'
#' @param ref A matrix
#' @param mix A matrix
#' @return A list including the ref and mix normalized and a ggplot
rnaseq_norm_ref_mix = function(ref, mix){

	mix_rn = rownames(mix)
	ref_rn = rownames(ref)
	
	if(max(ref, na.rm=T) < 50 | max(mix) < 50) stop("ARMET: Both objects have to be in count form, log tranformation detected")
	#print(head(mix))
	mix = rnaseq_norm(mix)
	#print(head(cbind(matrixStats:::rowMedians(mix), matrixStats:::rowMedians(ref, na.rm=T))))
	nf.obj = rnaseq_norm.calcNormFactor(cbind(matrixStats:::rowMedians(mix), matrixStats:::rowMedians(ref, na.rm=T)), 1)
	nf = nf.obj$nf
	my_df = nf.obj$df
	
	mix = mix / (sum(my_df[,1], na.rm=T) * nf[1]) *  sum(my_df[,1], na.rm = T)
	ref = ref / (sum(my_df[,2], na.rm=T) * nf[2]) *  sum(my_df[,1], na.rm = T)
	
	rownames(mix) = mix_rn 
	rownames(ref) = ref_rn 
	
	p = plot_densities(cbind(ref,mix))
	
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
quant_norm_to_RNAseq = function(target, obj){
	writeLines("Quantile norm to RNA seq..")
	
	# source("~/PhD/deconvolution/ARMET_dev/ARMET_FI.R")
	# res =  ARMET_FI(target)
	# target = res$df
	#if(!is.null(ncol(target))) if(ncol(target)>1) target = rnaseq_norm(target)
	
	do_log_mix =                     if(max(target, na.rm=T) > 50) T else F
	do_log_ref =                     if(max(obj) > 50) T else F
	
	if(do_log_mix) target =          log(target + 1)
	if(do_log_ref) obj =             log(obj + 1)
	
	target = as.matrix(target)
	obj = as.matrix(obj)
	
	rn = rownames(obj)
	cn = colnames(obj)
	obj = preprocessCore:::normalize.quantiles.use.target(obj, target=as.vector(target))
	rownames(obj) = rn
	colnames(obj) = cn
	
	if(do_log_ref) obj =             exp(obj) -1
	if(do_log_mix) target =          exp(target) -1
	
	p = plot_densities_double(target,obj)
	
	return(list (ref = target, mix = obj, plot=p ))
}

#' Normalize RNAseq to match array
#'
#' @param target A matrix
#' @param obj A matrix
#' @return A list including the ref and mix normalized and a ggplot
quant_norm_to_array = function(obj, target){
	writeLines("Quantile norm..")
	
	target = as.matrix(target)
	obj = as.matrix(obj)
	
	do_log_mix =                     if(max(target) > 50) T else F
	do_log_ref =                     if(max(obj) > 50) T else F
	
	if(do_log_mix) target =          log(target + 1)
	if(do_log_ref) obj =             log(obj + 1)
	
	rn = rownames(obj)
	cn = colnames(obj)
	obj = preprocessCore:::normalize.quantiles.use.target(obj, target=as.vector(target))
	rownames(obj) = rn
	colnames(obj) = cn
	
	if(do_log_ref) obj =             exp(obj) - 1
	if(do_log_mix) target =             exp(target) - 1
	obj[obj<0] = 0
	target[target<0] = 0
	
	p = plot_densities_double(target,obj)
	
	return(list (ref = target, mix = obj,  plot = p ))
}

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
		quant_norm_to_RNAseq(ref, mix)
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
prepare_input = function(ref, cell_types, markers){
	
	# Check if sd == 0 for any marker
	if(max(ref, na.rm=T)>0) ref = check_if_sd_zero(ref)
	
	e.obj = reshape:::melt(ref)
	e.obj$gene_num = match(e.obj$X1, markers) 
	e.obj$ct_num = match(e.obj$X2, cell_types)
	
	e_mu.obj = unique(e.obj[, c("gene_num", "ct_num")])
	
	e.obj$map_to_mu = match(paste(e.obj$gene_num, e.obj$ct_num), paste(e_mu.obj$gene_num, e_mu.obj$ct_num))
	
	return(list(e=e.obj, e_mu=e_mu.obj))
	
}

#' Averages the signatures of the same cell type
#'
#' @param df A matrix
#' @param do_log A bool
#' @param verbose A bool
#' @return A averaged matrix
get_mean_signature = function(df, do_log=F, verbose= T){
	if(nrow(df)==0) stop("ARMET: in get_mean_signature the data set has not records. Check bug upstream")

	if(max(df, na.rm=T)>50) do_log = T
	if(do_log) df = log(df+1)
	temp =                     data.frame(ct = colnames(df), t(df), row.names=NULL)
	df.mean =                  do.call("cbind", by(temp, temp$ct, function(df) {
		#print(unique(df[,1]))
		matrixStats:::colMedians(as.matrix(df[,-1]), na.rm=T)
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
