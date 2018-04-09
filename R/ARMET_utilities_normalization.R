#' Eliminate lowly tanscribed genes for normalization
#'
#' @param df A matrix
#' @param cpm_threshold A number
#' @param prop A number
#' @return A tibble filtered
rnaseq_norm.get_cpm = function(df,cpm_threshold = 0.5){
	
	cpm_threshold = 
		cpm_threshold / 
		(
			df %>%
				dplyr::group_by(sample) %>%
				dplyr::summarise(s = sum(value)) %>%
				dplyr::ungroup() %>%
				dplyr::summarise(m =median(s)) %>%
				dplyr::pull(m) /
				1e6
		)

	df %>% 
	{ if("ct"%in%names(.)) dplyr::select(-ct)	else .}	%>%
		dplyr::select(gene, sample, value) %>% 
		tidyr::spread(sample, value) %>% 
		tidyr::drop_na() %>%
		dplyr::do(
			dplyr::bind_cols(
				gene = (.)$gene,
				tibble::as_tibble( edgeR::cpm(	(.) %>% dplyr::select(-gene) ) )
			)
		) %>%
		tidyr::gather(sample, cpm, -gene) %>%
		dplyr::mutate(cpm_threshold = cpm_threshold)
	
}


#' Eliminate lowly tanscribed genes for normalization
#'
#' @param df A matrix
#' @param cpm_threshold A number
#' @param prop A number
#' @return A tibble filtered
rnaseq_norm.get_low_expressed = function(df,cpm_threshold = 0.5, prop = 3/4){

	cpm_threshold = 
		cpm_threshold / 
		(
			df %>%
				dplyr::group_by(sample) %>%
				dplyr::summarise(s = sum(value)) %>%
				dplyr::ungroup() %>%
				dplyr::summarise(m =median(s)) %>%
				dplyr::pull(m) /
				1e6
		)

	df %>% 
	{ if("ct"%in%names(.)) dplyr::select(-ct)	else .}	%>%
		dplyr::select(gene, sample, value) %>% 
		tidyr::spread(sample, value) %>% 
		tidyr::drop_na() %>%
		dplyr::do(
			(.) %>% 
				dplyr::filter(
					rowSums(
						edgeR::cpm(
							(.) %>% 
								dplyr::select(-gene)
						) > cpm_threshold
					) <
						floor(ncol( 
							(.) %>% dplyr::select(-gene)
						) * prop)
				)
		) %>% 
		droplevels() %>%
		dplyr::pull(gene) %>%
		as.character()
}

error_if_log_transformed = function(x){
	if(length(x$value)>0) if(max(x$value, na.rm=T)<50) 
		stop("ARMET: The input was log transformed in: check_if_sd_zero_and_correct")
}

#' Calculate the norm factor with calcNormFactor from limma
#'
#' @param df A matrix
#' @param reference A reference matrix
#' @param cpm_threshold A number
#' @param prop A number
#' @return A list including the filtered data frame and the normalization factors
rnaseq_norm.calcNormFactor = function(df, reference = NULL, cpm_threshold = 0.5, prop = 3/4, genes_to_keep = c()){
	
	error_if_log_transformed(df)
	
	# Get list of low transcribed genes
	gene_to_exclude = 
		rnaseq_norm.get_low_expressed(
			df %>% dplyr::filter(sample!="reference") ,
			cpm_threshold = cpm_threshold, 
			prop = prop
		)
	
	if(length(gene_to_exclude) == df %>% dplyr::distinct(gene) %>% nrow()) stop("ARMET: The gene expression matrix has been filtered completely for lowly expressed genes")
	
	# Keep genes that are forced 
	gene_to_exclude = gene_to_exclude[!gene_to_exclude %in% genes_to_keep]
	
	writeLines(sprintf("ARMET: %s genes excluded for normalization", length(gene_to_exclude)))
	
	df.filt =
		df %>%
		dplyr::filter(!gene %in% gene_to_exclude) %>%
		droplevels()
	
	list(
		gene_to_exclude = gene_to_exclude,
		nf = 
			tibble::tibble(
				sample = factor(levels(df.filt$sample)),
				nf = edgeR::calcNormFactors(
					df.filt %>% 
						tidyr::spread(sample, value) %>% 
						dplyr::select(-gene), 
					refColumn=reference, 
					method="TMM"
				),
			) %>%
			dplyr::left_join(
				df.filt %>%
					dplyr::group_by(sample) %>%
					dplyr::summarise(tot_filt = sum(value, na.rm = T)) %>%
					dplyr::mutate(sample = as.factor(as.character(sample))),
				by = "sample"
			)
	)
}

#' Normalize a RNA seq data set using rnaseq_norm.calcNormFactor
#'
#' @param df A matrix
#' @param reference A reference matrix
#' @param cpm_threshold A number
#' @param prop A number
#' @return A list including the filtered data frame and the normalization factors
rnaseq_norm = function(df, reference = NULL, cpm_threshold = 0.5, prop = 3/4, genes_to_keep = c()){
	#if(verbose) writeLines("Normalizing RNA-seq data with TMM")
	
	if(length(intersect( c("sample", "gene"), colnames(df))) < 2 ) stop("ARMET: input table is not correctly formatted as gene, sample")

	# Get norm factor object
	nf_obj = rnaseq_norm.calcNormFactor(df, reference, cpm_threshold, prop, genes_to_keep = genes_to_keep)
	
	# Calculate normalization factors
	nf = nf_obj$nf %>%
		dplyr::left_join(
			df %>% 
				dplyr::group_by(sample) %>% 
				dplyr::summarise(tot = sum(value, na.rm = T)) %>%
				dplyr::mutate(sample = as.factor(as.character(sample))),
			by = "sample"
		) %>%
		dplyr::mutate(multiplier = 1 / (tot_filt * nf) * max(tot) ) %>%
		{
			# I have correct the strange behaviour of edgeR of reference 
			# sample not being 1
			if("reference" %in% ( (.) %>% dplyr::pull(sample) ))
				(.) %>% 
				dplyr::mutate(
					multiplier = 
						multiplier / 
						(.) %>% 
						dplyr::filter(sample == "reference") %>% 
						dplyr::pull(multiplier)
				) 
			else
				(.)
		}
		
	df %>% 
		dplyr::mutate(sample = as.factor(as.character(sample))) %>%
		dplyr::left_join( nf, by = "sample") %>%
		dplyr::mutate(
			value_original = value, 
			value = value * multiplier
		) %>%
		dplyr::mutate(filtered_out = gene %in% nf_obj$gene_to_exclude) %>%
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
	rnaseq_norm(
		dplyr::bind_rows(
			target %>% 
				dplyr::group_by(gene) %>% 
				dplyr::summarise(value = median(value, na.rm=T)) %>% 
				dplyr::mutate(sample="reference"),
			obj %>% 
				dplyr::mutate(sample = as.character(sample))
		) %>% 
		# Put sample reference first for calc norm factor
		dplyr::mutate(
			sample = factor(
				sample,
				levels=c("reference", unique(sample[sample!="reference"]))
			)
		),
		reference = 1
	) %>%
	dplyr::filter(sample!="reference") %>%
	droplevels()

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
	if(!is_mix_microarray) rnaseq_norm_ref_mix(mix, ref)

	else {
		norm = quant_norm_to_target(mix, ref)
		list(mix=norm$obj, ref=norm$target)
	}
	
	
}
