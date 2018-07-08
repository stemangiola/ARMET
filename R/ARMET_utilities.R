
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
	library(jsonlite)
	tree = read_json("/wehisan/home/allstaff/m/mangiola.s/PhD/deconvolution/ARMET_dev/ARMET_TME_tree_RNAseq.json")
	save(tree, file="data/tree_json.rda")

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

error_if_log_transformed = function(x){
	if(length(x$value)>0) if(max(x$value, na.rm=T)<50)
		stop("ARMET: The input was log transformed in: check_if_sd_zero_and_correct")
}

#' Check the input for anomalies
check_input = function(mix, is_mix_microarray, my_design, cov_to_test, prior_sd, custom_ref, tree){

	# Check if mix is a tibble
	if(!tibble::is_tibble(mix)) stop("ARMET: The mixture must be a tibble")

	# Check if custom ref is tibble
	if(!is.null(custom_ref)){
		if(!tibble::is_tibble(custom_ref))
			stop("ARMET: The reference must be a tibble")
		if(!all(levels(ref_RNAseq_recursive_summary$ct) %in%colnames(custom_ref)))
			stop(
				sprinf(
					"ARMET: The columns of the provided reference must be %s",
					paste(levels(ref_RNAseq_recursive_summary$ct), collapse=" ")
				)
			)
	}

	# Check columns of design matrix
	if(!is.null(my_design) && !"sample" %in% colnames(my_design))
		stop("ARMET: The design matrix should have a column called 'sample' being a factor")

	# Check if cov to check in design
	if(!is.null(cov_to_test) && !cov_to_test%in%colnames(my_design))
		stop("ARMET: cov_to_test should be among the column names of the design matrix")

	# Check if microarray data
	if(is_mix_microarray) writeLines("ARMET: The input matrix is from microarray.")

	# Check if NA in mix
	if(
		mix %>%
		dplyr::select(-gene) %>%
		dplyr::select_if(function(.) any(is.na(.))) %>%
		ncol() > 0
	) stop("ARMET: NAs found in the query matrix")

	# Check if NA in mix genes
	if(
		mix %>%
		dplyr::select(gene) %>%
		dplyr::select_if(function(.) any(is.na(.))) %>%
		ncol() > 0
	) stop("ARMET: NAs found in the gene names")


	# Check if negative numbers in the matrix
	if(
		mix %>%
		dplyr::select(-gene) %>%
		dplyr::select_if(function(.) any((.)<0)) %>%
		ncol() > 0
	) stop("ARMET: negative numbers found in query matrix")


	# Check if design is tibble
	if(!is.null(cov_to_test) && !tibble::is_tibble(my_design)) stop("ARMET: The design matrix must be a tibble")

	# # This is how many conditions are in the study (e.g., treatment-vs-non-treatment)
	# if(!is.null(my_design)) {
	# 	writeLines("ARMET: The design matrix is :")
	# 	print(my_design)
	# }

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
			ct_to_consider = unlist(c(h, get_leave_label(get_node_from_name(node, h), recursive = F)))

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
				count = length(which(get_leave_label(get_node_from_name(tree, ct), label = "markers") %in%	ref$gene	)),
				val = "gene"
			)
		},
		tibble::as_tibble(table(ref$ct)) %>%
			dplyr::rename(ct=Var1, count = n) %>%
			dplyr::mutate(val = "sample")
	)

}

# plot_model_results = function(){
#
# 	proportions_4_plot = as.data.frame(proportions)
# 	proportions_4_plot$sample = rownames(proportions_4_plot)
# 	proportions_4_plot$cov = if(length(cov_to_test)>0) my_design[,"cov"] else 1
# 	proportions_4_plot = reshape::melt(proportions_4_plot, id.vars=c("sample", "cov"))
# 	proportions_4_plot$perc_1 = reshape::melt(proportions_2.5)$value
# 	proportions_4_plot$perc_2 = reshape::melt(proportions_97.5)$value
# 	proportions_4_plot$minus_sd = proportions_4_plot$value - reshape::melt(proportions_sd)$value
# 	proportions_4_plot$plus_sd = proportions_4_plot$value + reshape::melt(proportions_sd)$value
# 	p = ggplot2::ggplot( proportions_4_plot, ggplot2::aes(x=jitter(cov), y=value,fill=factor(variable)))+
# 		ggplot2::geom_boxplot(coef = 6) +
# 		ggplot2::geom_point(size=0.1) +
# 		ggplot2::geom_linerange(ggplot2::aes(ymin = perc_1, ymax = perc_2),  alpha=0.05) +
# 		ggplot2::geom_errorbar(ggplot2::aes(ymin = minus_sd, ymax = plus_sd),width=.2,  alpha=0.2) +
# 		ggplot2::facet_grid(~variable) +
# 		ggplot2::theme(
# 			axis.text.x= ggplot2::element_text(angle=90, vjust=0.4,hjust=1),
# 			panel.background = ggplot2::element_blank()
# 		)
#
# 	#if(save_report) ggplot2::ggsave(sprintf("%s/%s_proportions.pdf", output_dir, ct), useDingbats=FALSE,plot=p)
#
# 	if(length(cov_to_test)>0){
# 		#plot generate quantities
# 		gq <- as.matrix(fit, pars = c("beta_gen"))[1,]
# 		gq= parse_summary_vector_in_2D(gq)
# 		colnames(gq) = colnames(mix)
# 		rownames(gq) = colnames(model.in$x)
# 		gq = as.data.frame(t(gq))
# 		gq$cov = my_design[,"cov"]
# 		gq$sample = rownames(gq)
# 		gq = reshape::melt(gq, id.vars=c("sample", "cov"))
# 		p = ggplot2::ggplot( gq, ggplot2::aes(x=factor(cov), y=value,fill=factor(variable)))+
# 			ggplot2::geom_boxplot()+
# 			ggplot2::geom_jitter(position=position_dodge(width=0.75), ggplot2::aes(group=variable)) +
# 			ggplot2::facet_grid(~variable) +
# 			ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, vjust=0.4,hjust=1))
#
# 		#if(save_report) ggplot2::ggsave(sprintf("%s/%s_proportions_posterior.pdf", output_dir, ct), useDingbats=FALSE,plot=p)
# 	}
# 	# Calculate predicted values
# 	#pred.df=as.matrix(ref[markers,])%*%t(as.matrix(res.df))
#
#
# 	# about means
# 	df_4_plot = data.frame(
# 		merge(
# 			reshape::melt(mix+1),
# 			reshape::melt(t(proportions %*% t(ref.mean)+1)),
# 			by=c("X1", "X2")
# 		)
# 	)
# 	colnames(df_4_plot)[3:4] = c("value", "predicted")
# 	df_4_plot$ct = sapply(df_4_plot$X1, function(x1) node_from_gene(my_tree, x1)$name )
# 	df_4_plot$bg = apply(df_4_plot, 1, function(mr) 	y_hat_background[mr[2], mr[1]] )
# 	df_4_plot$bg_prop = df_4_plot$bg /df_4_plot$value
# 	df_4_plot$bg_prop[df_4_plot$bg_prop >10] = 10
#
#
# 	df_4_plot = df_4_plot[df_4_plot$X2%in%unique(df_4_plot$X2)[1:20], ]
# 	#df_4_plot = df_4_plot[grep("CD8", df_4_plot$X2), ]
#
# 	p = ggplot2::ggplot(df_4_plot, ggplot2::aes(predicted, value, label = X1, size=bg_prop,color = factor(ct))) +
# 		ggplot2::geom_point(alpha=0.5) +
# 		ggplot2::scale_size(range = c(0, 5)) +
# 		ggplot2::expand_limits(x = 0.1, y = 0.1) +
# 		ggplot2::geom_abline(intercept=0, slope=1) +
# 		ggplot2::scale_x_log10() +
# 		ggplot2::scale_y_log10() +
# 		ggrepel::geom_text_repel(
# 			ggplot2::aes(color = ct),
# 			size = 2,
# 			segment.alpha= 0.2) +
# 		ggplot2::facet_wrap( ~ X2, nrow=5) +
# 		ggplot2::coord_fixed() +
# 		ggplot2::theme(
# 			panel.background = ggplot2::element_blank(),
# 			legend.text=ggplot2::element_text(size=30)
# 		)
#
# 	#if(save_report) ggplot2::ggsave(sprintf("%s/%s_predicted_values.pdf", output_dir, ct), useDingbats=FALSE, width = 80, height = 80, units = "cm", plot=p)
#
#
# 	if(!is.null(observed_prop) & 0) {
# 		writeLines("ARMET: start plotting expected gene vallues")
# 		link = tibble::as_tibble(do.call("rbind", lapply(unique(colnames(observed_prop)), function(op) {
# 			ct = get_hierarchy(my_tree, ct=op)
# 			c(op, ct[ct%in%colnames(ref.mean)])
# 		})))
# 		if(ncol(link)>=2){
# 			colnames(link) = c("ct", "ct_new")
# 			observed_prop = tibble::as_tibble(reshape::melt(observed_prop))
# 			colnames(observed_prop) = c("sample", "ct", "value")
#
# 			observed_prop =
# 				dplyr::left_join(observed_prop, link, by="ct") %>%
# 				dplyr::select(-ct) %>%
# 				dplyr::group_by(sample, ct_new) %>%
# 				dplyr::summarise(value = sum(value, na.rm=TRUE)) %>%
# 				tidyr::spread(ct_new, value)
#
# 			for(r in colnames(ref.mean)[!colnames(ref.mean)%in%colnames(observed_prop)]) {
# 				my_df = matrix(rep(0, nrow(observed_prop)))
# 				colnames(my_df) = r
# 				observed_prop = observed_prop %>% dplyr::bind_cols(as_tibble(my_df))
# 			}
#
# 			observed_prop = as.data.frame(observed_prop)
# 			rownames(observed_prop) = observed_prop[,1]
# 			observed_prop = observed_prop[,-1]
# 			observed_prop = as.matrix(observed_prop)
# 			observed_prop = observed_prop[,colnames(ref.mean)]
#
# 			print(proportions)
#
# 			df_4_plot = data.frame(merge(reshape::melt(mix+1), reshape::melt(t(observed_prop %*% t(ref.mean)+1)), by=c("X1", "X2")))
# 			colnames(df_4_plot)[3:4] = c("value", "predicted")
# 			df_4_plot$ct = sapply(df_4_plot$X1, function(x1) node_from_gene(my_tree, x1)$name )
# 			df_4_plot$bg = apply(df_4_plot, 1, function(mr) 	y_hat_background[mr[2], mr[1]] )
# 			df_4_plot$bg_prop = df_4_plot$bg /df_4_plot$value
# 			df_4_plot$bg_prop[df_4_plot$bg_prop >10] = 10
#
# 			p = ggplot2::ggplot( df_4_plot[df_4_plot$X2%in%unique(df_4_plot$X2)[1:20], ], ggplot2::aes(predicted, value, label = X1, size=bg_prop)) +
# 				ggplot2::geom_point(alpha=0.5) +
# 				ggplot2::scale_size(range = c(0, 5)) +
# 				ggplot2::expand_limits(x = 0.1, y = 0.1) +
# 				ggplot2::geom_abline(intercept=0, slope=1) +
# 				ggplot2::scale_x_log10() +
# 				ggplot2::scale_y_log10() +
# 				ggrepel::geom_text_repel(
# 					aes(color = ct),
# 					size = 2,
# 					segment.alpha= 0.2) +
# 				ggplot2::facet_wrap( ~ X2, nrow=5) +
# 				ggplot2::coord_fixed() +
# 				ggplot2::theme(
# 					panel.background = ggplot2::element_blank(),
# 					legend.text=ggplot2::element_text(size=30)
# 				)
#
# 			#if(save_report) ggplot2::ggsave(sprintf("%s/%s_predicted_values_observed_prop.pdf", output_dir, ct), useDingbats=FALSE, width = 80, height = 80, units = "cm", plot=p)
# 		}
#
# 	}
#
# }


#' Parse the annotated tree
#' @rdname ARMET_getFit
#'
#' Prints a report of the hipothesis testing
#'
#' @param ARMET-tc object
#'
#' @return a print out of the hierarchy
#'
#' @examples
#'  ARMET_getFit(ARMET_tc_result)
#' @export
ARMET_getFit = function(obj){
	print(
		obj$stats,
		"Estimate" ,
		"Direction",
		"CI.low" ,
		"CI.high" ,
		"Sig" ,
		"Driver"
	)
}

ref_to_summary_ref = function(tree, ref){

	# cl <- multidplyr::create_cluster(parallel::detectCores()-2)
	# multidplyr::set_default_cluster(cl)

	get_distributions_from_tree(
		get_gene_distributons_recursive(
			tree,
			# Get full matrix expression
			ref %>%
				#Get summary statistics in parallel
				#multidplyr::partition(ct, gene) %>%
				dplyr::group_by(ct, gene) %>%
				dplyr::summarise(
					log_mean = median(log(value+1), na.rm = T),
					log_sd =   mad(log(value+1), na.rm = T)
				) %>%
				#dplyr::collect() %>%
				dplyr::ungroup()
		)
	) %>%
		dplyr::mutate(value = exp(log_mean))

	# parallel::stopCluster(cl)
}

# Add hypothesis testing
parse_fit_for_quantiles = function(fit, param = "beta", q, label, my_design, names_groups){

	as.data.frame(rstan::summary(fit, param, probs = q)$summary) %>%
		tibble::as_tibble(rownames="par") %>%
		tidyr::separate(par, c("dummy", "sample_num", "ct_num", "dummy2"), "\\[|,|\\]") %>%
		dplyr::mutate(ct_num = as.integer(ct_num), sample_num = as.integer(sample_num)) %>%
		dplyr::select(-dummy, -dummy2) %>%
		dplyr::filter(ct_num > max( ct_num ) -length(names_groups)) %>%
		dplyr::select(1:2,ifelse(q == 0.5, "mean", ncol(.) - 2)) %>%
		tidyr::spread(2,3) %>%
		dplyr::select(-sample_num) %>%
		stats::setNames(names_groups) %>%
		dplyr::mutate(sample = levels(my_design$sample)) %>%
		tidyr::gather(ct, !!label, -sample) %>%
		dplyr::mutate_if(is.character, as.factor)

}

any_column_double = function(x) is.numeric(x) & !all(x %in% c(0:1))

annotate_posterior_sample_ct = function(tbl, ct_names, my_design) {
	tbl %>%
	left_join(
		tibble( ct_idx = 1:length(ct_names), ct=ct_names ),
		by="ct_idx"
	) %>%
	left_join(
		bind_cols( sample_idx = 1:nrow(my_design), my_design %>% select(-`(Intercept)`) %>% arrange(sample) ),
		by="sample_idx"
	)
}

