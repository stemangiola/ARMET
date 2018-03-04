#' # Functions for pvalue calc
#' test_two_distrib = function(mu_1, sd_1, mu_2, sd_2) 2*pnorm( -abs(mu_1-mu_2)/(sum_euclidean(c(sd_1, sd_2))))
#' test_point_distrib = function(point, mu, sd) 2*pnorm(point, mu, sd, lower.tail = mu > point)
#' sum_euclidean = function(x) sqrt(sum(x^2))
#' 
#' get_row_col_from_unique_index = function(mat, i){
#' 	my_row = floor(i/ncol(mat))
#' 	if(i-(my_row*ncol(mat)) == 0) {
#' 		my_col =  ncol(mat)
#' 	} else {
#' 		my_col =  i - my_row*ncol(mat)
#' 		my_row = my_row + 1
#' 	}
#' 	
#' 	list(ro = my_row, co = my_col)
#' }
#' 
#' get_closest_pair = function(cap){
#' 	
#' 	if(!all(c("sample", "value") %in% colnames(cap))) stop("ARMET: Colnames_do_not_match")
#' 	
#' 	my_test = function(pair, cap, check_background ){
#' 		
#' 		# Eliminate combination that are not with background, if present
#' 		if(check_background && !"background" %in% pair) return(NA)
#' 		
#' 		test_two_distrib(
#' 			mean((cap %>% dplyr::filter(sample==pair[1]))$value),
#' 			sd((cap %>% dplyr::filter(sample==pair[1]))$value),
#' 			mean((cap %>% dplyr::filter(sample==pair[2]))$value),
#' 			sd((cap %>% dplyr::filter(sample==pair[2]))$value)
#' 		)
#' 	} 
#' 	
#' 	# Get closest
#' 	my_samples = cap %>% dplyr::distinct(sample) %>% dplyr::pull(sample)
#' 	my_pvalues = combn(x = my_samples, m = 2, simplify = TRUE, FUN= my_test,  cap, "background"%in%my_samples)
#' 	as.data.frame(t(combn(my_samples, 2))) %>% 
#' 		stats::setNames(c("s1", "s2")) %>%
#' 		tibble::as_tibble() %>%
#' 		dplyr::bind_cols(p = my_pvalues) %>%
#' 		dplyr::filter(p == max(p, na.rm = T)) %>%
#' 		mutate_if(is.factor, as.character)%>%
#' 		tidyr::gather(value, sample, -p) %>%
#' 		select(-value)
#' }
#' 
#' get_pvalues_from_distributons = function(coef_ang_posterior){
#' 	
#' 	cap = 
#' 		as.data.frame(coef_ang_posterior) %>% 
#' 		tibble::as_tibble() %>%
#' 		tidyr::gather(sample, value)
#' 	background = cap[0,]
#' 	significance = tibble::tibble(p = numeric(), sample = character())
#' 	
#' 	
#' 	# Iterate untill all parameters have been tested
#' 	while(nrow(cap)>0){
#' 		closest_pair = get_closest_pair(rbind( cap, background) )
#' 		significance = significance %>% dplyr::bind_rows( closest_pair %>% dplyr::filter(sample != "background") )
#' 		
#' 		background = background %>%
#' 			bind_rows(
#' 				cap %>%	
#' 					dplyr::filter(
#' 						# Samples that are in the cluster
#' 						sample %in% (
#' 							closest_pair %>% 
#' 								dplyr::filter(p > 0.05) %>% 
#' 								dplyr:: pull(sample)
#' 						)
#' 					)
#' 			) %>%
#' 			# If the first cluster is not even a cluster 
#' 			# create a sunthetic one
#' 			{
#' 				if(nrow(background) == 0)
#' 					cap %>%
#' 					dplyr::filter( sample %in% ( closest_pair %>% dplyr:: pull(sample) ))  %>% 
#' 					dplyr::mutate(
#' 						value = rnorm(
#' 							n(), 
#' 							mean(value),
#' 							# SD is the mean of the SDs of the closest parameters
#' 							(.) %>% 
#' 								dplyr::group_by(sample) %>% 
#' 								dplyr::summarise(s = sd(value)) %>% 
#' 								dplyr::pull(s) %>% 
#' 								mean()
#' 						)
#' 					)
#' 				else (.)
#' 			} %>%
#' 			# Change name to background
#' 			mutate(sample = "background")
#' 		
#' 		cap = cap %>% dplyr::filter(!sample %in% significance$sample)
#' 		
#' 	}
#' 	
#' 	list(significance = significance, center_bg = mean(background$value))
#' }
#' 
#' 
#' #' hypothesis test for the covariate of choice
#' #' @rdname dirReg_test
#' #'
#' #' Prints a report of the hipothesis testing
#' #'
#' #' @param fit stan fit object 
#' #' @param my_design design matrix
#' #' @param cov_to_test character variable indicating the column name of the design matrix to test
#' #' @param which_cov_to_test for internal usage
#' #'
#' #' @return a vector including
#' #'     pvalues of the intercepts
#' #'     pvalues of the angular coefficients
#' #'     sign of the angular coefficient
#' #'
#' #' @examples
#' #'  dirReg_test(fit, my_design, cov_to_test)
#' #' @export
#' dirReg_test = function(fit, my_design, cov_to_test = NULL, which_cov_to_test = 2, names_groups = NULL){
#' 	
#' 	logit_adj <- function(v, t=0.5) -log(t*(v-1) / ((t-1)*v));
#' 	
#' 	# Decide which covariate to check
#' 	if(!is.null(cov_to_test)) which_cov_to_test = which(colnames( my_design %>% dplyr::select(-sample) ) == cov_to_test )
#' 	
#' 	# Get the posterior distributions
#' 	alpha_posterior = as.matrix(as.data.frame(rstan:::extract( fit, "alpha")))
#' 	coef_ang_posterior = as.matrix(alpha_posterior[,grep(sprintf("alpha.%s", which_cov_to_test), colnames(alpha_posterior), fixed=T)])
#' 	interc_posterior = as.matrix(alpha_posterior[,grep("alpha.1", colnames(alpha_posterior), fixed=T)])
#' 	
#' 	
#' 	
#' 	K = ncol(interc_posterior)
#' 	
#' 	obj = get_pvalues_from_distributons(coef_ang_posterior)
#' 	p = obj$significance
#' 	m = obj$center_bg
#' 	
#' 	coef_ang_posterior_adj = 
#' 		data.frame(logit_adj(coef_ang_posterior, m)) %>%
#' 		tibble::as_tibble() %>%
#' 		stats::setNames(names_groups) %>%
#' 		tidyr::gather(ct, value)
#' 	
#' 
#' 	
#' 	stats = do.call("rbind", lapply(1:K, function(i){
#' 		
#' 		# Angular coefficient
#' 		mcap = mean(coef_ang_posterior[,i])
#' 		scap = sd(coef_ang_posterior[,i])
#' 		ecap = scap/sqrt(nrow(coef_ang_posterior))
#' 		pcap = p %>% dplyr::filter(sample == colnames(coef_ang_posterior)[i]) %>% dplyr::pull(p)
#' 		mcap_human_readable = logit_adj(mcap, m)
#' 		pcap_human_readable = if(pcap<0.01 & pcap>0) formatC(pcap, format = "e", digits = 2) else round(pcap, 2)
#' 		
#' 		# ntercept
#' 		mip = mean(interc_posterior[,i])
#' 		sip = sd(interc_posterior[,i])
#' 		eip = sip/sqrt(nrow(interc_posterior))
#' 		pip = test_point_distrib(1/K, mip, sip)  
#' 		mip_human_readable = mip
#' 		pip_human_readable = if(pip<0.01& pip>0) formatC(pip, format = "e", digits = 2) else round(pip, 2)
#' 		
#' 		tibble::tibble(
#' 			m = m,
#' 			mcap = mcap,
#' 			scap = scap,
#' 			ecap = ecap,
#' 			mip = mip,
#' 			sip = sip,
#' 			eip = eip,
#' 			pcap = pcap,
#' 			pip = pip,
#' 			mip_human_readable = mip_human_readable,
#' 			mcap_human_readable = mcap_human_readable,
#' 			pcap_human_readable = pcap_human_readable,
#' 			pip_human_readable = pip_human_readable,
#' 			symbol_pcap = if(pcap<0.001) "***" else if(pcap<0.01) "**" else if(pcap<0.05) "*" else if(pcap<0.1) "." else "",
#' 			symbol_pip = if(pip<0.001) "***" else if(pip<0.01) "**" else if(pip<0.05) "*" else if(pip<0.1) "." else "",
#' 			direction = if(mcap_human_readable>0) "+" else if(mcap_human_readable<0) "-" else "0"
#' 		)
#' 		
#' 	})) %>%
#' 		dplyr::mutate(ct = names_groups ) %>%
#' 		dplyr::mutate_if(is.character, as.factor)
#' 	
#' 	list(
#' 		stats=stats, 
#' 		estimate_prop_with_uncertanties = estimate_prop_with_uncertanties,
#' 		coef_ang_posterior_adj = coef_ang_posterior_adj
#' 	)
#' }
#' 
#' 
#' 
#' parse_extract_2D <- function(fit, param, fun)	apply(rstan::extract(fit, param)[[1]], c(2,3), fun )
#' 
#' #' @rdname beta_reg_hierarchical
#' #' @export
#' beta_reg_hierarchical = function(fit, my_design){
#' 	
#' 	writeLines("ARMET: Starting inference beta reg..")
#' 	
#' 	beta_posterior_mu = parse_extract_2D(fit, "beta", mean)
#' 	beta_posterior_sd = parse_extract_2D(fit, "beta", sd)
#' 	beta_posterior = aperm(rstan::extract(fit, "beta")[[1]], c(2,3,1))
#' 	beta_posterior = beta_posterior[,,1:100]
#' 	
#' 	rstan:::sampling( stanmodels$beta_reg_hierarchical, 
#' 										#stan(file="~/PhD/simplexRegression/src/stan_files/beta_reg_hierarchical.stan",
#' 										data=list(
#' 											I = dim(beta_posterior)[3],
#' 											K=ncol(beta_posterior_mu),
#' 											N=nrow(beta_posterior_mu),
#' 											X=my_design,
#' 											R=ncol(my_design),
#' 											beta_posterior = beta_posterior,
#' 											beta_mu=beta_posterior_mu,
#' 											beta_sd = beta_posterior_sd
#' 										),
#' 										cores=4,
#' 										iter = 1000,
#' 										refresh = 0
#' 	)
#' }
#' 
#' #' hypothesis test for the covariate of choice
#' #' @rdname betaReg_test
#' #'
#' #' Prints a report of the hipothesis testing
#' #'
#' #' @param fit stan fit object 
#' #' @param my_design design matrix
#' #' @param cov_to_test character variable indicating the column name of the design matrix to test
#' #' @param which_cov_to_test for internal usage
#' #'
#' #' @return a vector including
#' #'     pvalues of the intercepts
#' #'     pvalues of the angular coefficients
#' #'     sign of the angular coefficient
#' #'
#' #' @examples
#' #'  betaReg_test(fit, my_design, cov_to_test)
#' #' @export
#' betaReg_test = function(fit, my_design, cov_to_test = NULL, which_cov_to_test = 2, names_groups = NULL){
#' 	
#' 	# Decide which covariate to check
#' 	if(!is.null(cov_to_test)) which_cov_to_test = which(colnames( my_design %>% dplyr::select(-sample) ) == cov_to_test )
#' 	
#' 	# Get the posterior distributions
#' 	alpha_posterior = as.data.frame(rstan:::extract( fit, "alpha2"))
#' 	coef_ang_posterior = as.matrix(alpha_posterior[,grep(sprintf("alpha2.%s", which_cov_to_test), colnames(alpha_posterior), fixed=T)])
#' 	interc_posterior = as.matrix(alpha_posterior[,grep("alpha2.1", colnames(alpha_posterior), fixed=T)])
#' 	
#' 	K = ncol(interc_posterior)
#' 	m = 0
#' 	
#' 	plot = plot_densities(coef_ang_posterior , do_log = F, color="0",  fill = 1:K) +
#' 		ggtitle("Effectual trends - posterior distribution of angular coeff. - Beta regression")
#' 	
#' 	stats = 	do.call("rbind", lapply(1:K, function(i){
#' 		
#' 		mcap = mean(coef_ang_posterior[,i])
#' 		scap = sd(coef_ang_posterior[,i])
#' 		ecap = scap/sqrt(nrow(alpha_posterior))
#' 		mip = mean(interc_posterior[,i])
#' 		sip = sd(interc_posterior[,i])
#' 		eip = sip/sqrt(nrow(alpha_posterior))
#' 		pcap = min( 2*(1-pnorm(m, mcap, scap)), 2*(1-pnorm(m, mcap, scap, lower.tail = F)))
#' 		pip = min ( 2*(1-pnorm(1/K, mip, sip)), 2*(1-pnorm(1/K, mip, sip, lower.tail = F)))
#' 		mip_human_readable = inv_logit(mip)
#' 		mcap_human_readable = mcap
#' 		pcap_human_readable = if(pcap<0.01 & pcap>0) formatC(pcap, format = "e", digits = 2) else round(pcap, 2)
#' 		pip_human_readable = if(pip<0.01& pip>0) formatC(pip, format = "e", digits = 2) else round(pip, 2)
#' 		
#' 		tibble(
#' 			m = m,
#' 			mcap = mcap,
#' 			scap = scap,
#' 			ecap = ecap,
#' 			mip = mip,
#' 			sip = sip,
#' 			eip = eip,
#' 			pcap = pcap,
#' 			pip = pip,
#' 			mip_human_readable = mip_human_readable,
#' 			mcap_human_readable = mcap_human_readable,
#' 			pcap_human_readable = pcap_human_readable,
#' 			pip_human_readable = pip_human_readable,
#' 			symbol_pcap = if(pcap<0.001) "***" else if(pcap<0.01) "**" else if(pcap<0.05) "*" else if(pcap<0.1) "." else "",
#' 			symbol_pip = if(pip<0.001) "***" else if(pip<0.01) "**" else if(pip<0.05) "*" else if(pip<0.1) "." else "",
#' 			direction = if(mcap_human_readable>0) "+" else if(mcap_human_readable<0) "-" else "0"
#' 		)
#' 		
#' 	}))
#' 	stats$ct = names_groups
#' 	
#' 	stats = stats %>% dplyr::mutate_if(is.character, as.factor)
#' 	
#' 	list(stats=stats, plot=plot)
#' 	
#' }
