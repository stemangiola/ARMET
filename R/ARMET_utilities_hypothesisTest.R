
# parallel -j 10 'Rscript ~/PhD/test_phi.R {} | grep "^__" >> phi_test_10_6.csv' ::: 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54 56 58 60 62 64 66 68 70 72 74 76 78 80 82 84 86 88 90 92 94 96 98 100 102 104 106 108 110 112 114 116 118 120 122 124 126 128 130 132 134 136 138 140 142 144 146 148 150 152 154 156 158 160 162 164 166 168 170 172 174 176 178 180 182 184 186 188 190 192 194 196 198 200
# args = commandArgs(trailingOnly=TRUE)
# n = as.numeric(args[1])
# 
# library(devtools)
# library(simplexRegression)
# #devtools::load_all("~/PhD/simplexRegression")
# 
# n_samples = 300
# set.seed(123)
# my_design = model.matrix(~cov , data = data.frame(cov=runif(n_samples)))
# #my_design = model.matrix(~cov , data = data.frame(cov=sample(0:1, size=n_samples, replace=T)))
# 
# obj = simulate_proportions(
# 	n ,
# 	c(0.35, 0.1, 0.1, 0.1 , 0.1, 0.1) ,
# 	c(0, 8, 0, 0, 0, 0),
# 	my_design,
# 	n_samples = n_samples
# )
# 
# fit = simplexRegression(obj$prop, my_design, cov_to_test="cov")
# writeLines(sprintf("__ %s %s", n, paste(summary(fit, c("phi", "phi2"))$summary[,1], collapse=" ")))

parse_extract_2D <- function(fit, param, fun){
	apply(extract(fit, param)[[1]], c(2,3), fun )
}


logsumexp <- function (x) {
	y = max(x)
	y + log(sum(exp(x - y)))
}

softmax <- function (x) {
	exp(x - logsumexp(x))
}

logit <- function(v) log(v / (1 - v));

logit_adj <- function(v, t=0.5) -log(t*(v-1) / ((t-1)*v));

inv_logit <- function (u) 1 / (1 + exp(-u));


# Functions for pvalue calc
sum_euclidean = function(x) sqrt(sum(x^2))
test_two_distrib = function(mu_1, sd_1, mu_2, sd_2) 2*pnorm( -abs(mu_1-mu_2)/(sum_euclidean(c(sd_1, sd_2))))
test_point_distrib = function(point, mu, sd) 2*pnorm(point, mu, sd, lower.tail = mu > point)
test_point_distrib_CI = function(point, distrib) 2 * ecdf(distrib)(point) %>% { if((.)>0.5) 1-(.) else (.) } 

get_row_col_from_unique_index = function(mat, i){
	my_row = floor(i/ncol(mat))
	if(i-(my_row*ncol(mat)) == 0) {
		my_col =  ncol(mat)
	} else {
		my_col =  i - my_row*ncol(mat)
		my_row = my_row + 1
	}
	
	list(ro = my_row, co = my_col)
}

get_closest_pair = function(cap){
	
	if(!all(c("sample", "value") %in% colnames(cap))) stop("ARMET: Colnames_do_not_match")
	
	my_test = function(pair, cap, check_background ){
		
		# Eliminate combination that are not with background, if present
		if(check_background && !"background" %in% pair) return(NA)
		
		test_two_distrib(
			mean((cap %>% dplyr::filter(sample==pair[1]))$value),
			sd((cap %>% dplyr::filter(sample==pair[1]))$value),
			mean((cap %>% dplyr::filter(sample==pair[2]))$value),
			sd((cap %>% dplyr::filter(sample==pair[2]))$value)
		)
	} 
	
	# Get closest
	my_samples = cap %>% dplyr::distinct(sample) %>% dplyr::pull(sample)
	my_pvalues = combn(x = my_samples, m = 2, simplify = TRUE, FUN= my_test,  cap, "background"%in%my_samples)
	as.data.frame(t(combn(my_samples, 2))) %>% 
		stats::setNames(c("s1", "s2")) %>%
		tibble::as_tibble() %>%
		dplyr::bind_cols(p = my_pvalues) %>%
		dplyr::filter(p == max(p, na.rm = T)) %>%
		dplyr::mutate_if(is.factor, as.character)%>%
		tidyr::gather(value, sample, -p) %>%
		dplyr::select(-value)
}

get_pvalues_from_distributons = function(coef_ang_posterior){
	
	cap = 
		as.data.frame(coef_ang_posterior) %>% 
		tibble::as_tibble() %>%
		tidyr::gather(sample, value)
	background = cap[0,]
	significance = tibble::tibble(p = numeric(), sample = character())
	
	
	# Iterate untill all parameters have been tested
	while(nrow(cap)>0){
		closest_pair = get_closest_pair(rbind( cap, background) )
		significance = significance %>% dplyr::bind_rows( closest_pair %>% dplyr::filter(sample != "background") )
		
		background = background %>%
			dplyr::bind_rows(
				cap %>%	
					dplyr::filter(
						# Samples that are in the cluster
						sample %in% (
							closest_pair %>% 
								dplyr::filter(p > 0.05) %>% 
								dplyr:: pull(sample)
						)
					)
			) %>%
			# If the first cluster is not even a cluster 
			# create a sunthetic one
			{
				if(nrow(background) == 0)
					cap %>%
					dplyr::filter( sample %in% ( closest_pair %>% dplyr:: pull(sample) ))  %>% 
					dplyr::mutate(
						value = rnorm(
							n(), 
							mean(value),
							# SD is the mean of the SDs of the closest parameters
							(.) %>% 
								dplyr::group_by(sample) %>% 
								dplyr::summarise(s = sd(value)) %>% 
								dplyr::pull(s) %>% 
								mean()
						)
					)
				else (.)
			} %>%
			# Change name to background
			dplyr::mutate(sample = "background")
		
		cap = cap %>% dplyr::filter(!sample %in% significance$sample)
		
	}
	
	list(significance = significance, center_bg = mean(background$value))
}


#' hypothesis test for the covariate of choice
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
	
	obj = get_pvalues_from_distributons(coef_ang_posterior)
	p = obj$significance
	m = obj$center_bg
	
	coef_ang_posterior_adj = 
		data.frame(logit_adj(coef_ang_posterior, m)) %>%
		tibble::as_tibble() %>%
		stats::setNames(names_groups) %>%
		tidyr::gather(ct, value)
	
	stats = do.call("rbind", lapply(1:K, function(i){
		
		# Angular coefficient
		mcap = mean(coef_ang_posterior[,i])
		scap = sd(coef_ang_posterior[,i])
		ecap = scap/1
		#pcap = p %>% dplyr::filter(sample == colnames(coef_ang_posterior)[i]) %>% dplyr::pull(p)
		pcap = test_point_distrib_CI(m, coef_ang_posterior[,i])
		mcap_human_readable = logit_adj(mcap, m)
		pcap_human_readable = if(pcap<0.01 & pcap>0) formatC(pcap, format = "e", digits = 2) else round(pcap, 2)
		
		# ntercept
		mip = mean(interc_posterior[,i])
		sip = sd(interc_posterior[,i])
		eip = sip/1
		#pip = test_point_distrib(1/K, mip, sip)  
		pip = test_point_distrib_CI(1/K, interc_posterior[,i])
		mip_human_readable = mip
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
		
	})) %>%
		dplyr::mutate(ct = as.factor(names_groups )) 
	
	list(
		stats=stats, 
		coef_ang_posterior = coef_ang_posterior,
		coef_ang_posterior_adj = coef_ang_posterior_adj
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

	# Decide which covariate to check
	if(!is.null(cov_to_test)) which_cov_to_test = which(colnames(my_design)==cov_to_test)
	
	# Get the posterior distributions
	alpha_posterior = as.data.frame(rstan:::extract( fit, "alpha2"))
	coef_ang_posterior = as.matrix(alpha_posterior[,grep(sprintf("alpha2.%s", which_cov_to_test), colnames(alpha_posterior), fixed=T)])
	interc_posterior = as.matrix(alpha_posterior[,grep("alpha2.1", colnames(alpha_posterior), fixed=T)])
	
	K = ncol(interc_posterior)
	m = 0

	stats = 	do.call("rbind", lapply(1:K, function(i){
		
		mcap = mean(coef_ang_posterior[,i])
		scap = sd(coef_ang_posterior[,i])
		ecap = scap/sqrt(nrow(alpha_posterior))
		mip = mean(interc_posterior[,i])
		sip = sd(interc_posterior[,i])
		eip = sip/sqrt(nrow(alpha_posterior))
		#pcap = min( 2*(1-pnorm(m, mcap, scap)), 2*(1-pnorm(m, mcap, scap, lower.tail = F)))
		pcap = test_point_distrib_CI(m, coef_ang_posterior[,i])
		#pip = min ( 2*(1-pnorm(1/K, mip, sip)), 2*(1-pnorm(1/K, mip, sip, lower.tail = F)))
		pip = test_point_distrib_CI(1/K, interc_posterior[,i])
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
	stats$name = names_groups
	
	list(stats=stats, coef_ang_posterior = coef_ang_posterior)
	
}


#' @rdname beta_reg_hierarchical
#' @export
beta_reg_hierarchical = function(fit, my_design){
	
	writeLines("Starting inference beta reg..")
	
	beta_posterior_mu = parse_extract_2D(fit, "beta", mean)
	beta_posterior_sd = parse_extract_2D(fit, "beta", sd)
	beta_posterior = aperm(rstan:::extract(fit, "beta")[[1]], c(2,3,1))
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
	
#' Wrapper for the model and hipothesis testing of simplex regression
#'
#' @param proportions simplex matrix (rows must sum to 1)
#' @param my_design design matrix
#' @param cov_to_test character variable indicating the column name of the design matrix to test
#'
#' @return a stan fit object
#'
#' @examples
#'   simplexRegression(prop, my_design, "cov")
#'  
#' @rdname simplexRegression
#' @export
simplexRegression = function(proportions, my_design, cov_to_test = NULL, verbose = T){

 	writeLines("Starting inference second refined model..")

	# Check cov to test
	if(!is.null(cov_to_test) && !cov_to_test %in% colnames(my_design)) stop("cov_to_test must be in your design")
	
	# Correct for 0 prop ##############################
	###################################################

	#https://www.rdocumentation.org/packages/DirichletReg/versions/0.3-0/topics/DR_data
	proportions = ( proportions*(nrow(proportions)-1) + (1/ncol(proportions)) ) / nrow(proportions)
	# proportions[proportions==0] = min(proportions[proportions>0])
	# proportions = proportions/apply(proportions,1,sum)
	
	###################################################
	###################################################
	
	# Convert to data frame
	proportions = as.data.frame(proportions)
	my_design = as.data.frame(my_design) %>% tibble::as_tibble()
	
	if(ncol(my_design)>1) for(i in 2:ncol(my_design)) my_design[,i] = my_design[,i]/max(my_design[,i])
	
	fit = rstan::sampling( stanmodels$simplexRegression,
	  data=list(
	    K=ncol(proportions),
	    N=nrow(my_design),
	    X=my_design,
	    R=ncol(my_design),
	    beta=proportions
	  ),
	  cores=4,
	  refresh = 0
	)

	#traceplot(fit, pars="alpha")

	dt = dirReg_test(fit, my_design, cov_to_test, names_groups = colnames(proportions))
	bt = betaReg_test(fit, my_design, cov_to_test, names_groups = colnames(proportions))
	
	if(verbose){
		coef.mat = do.call("rbind", lapply(1:nrow(dt$stats), function(i){

			d = dt$stats[i,]
			b = bt$stats[i,]
			
			arr_names = c("name", "Estimate", "Std. Error", "z value", "Pr(>|z|)", "sig. annotation")
			arr_1 = c( "(Intercept)", round(d$mip, 2) , round(d$eip, 2), round((d$m - d$mip) / d$sip, 2), d$pip_human_readable, d$symbol_pip)
			arr_2 = c(paste(cov_to_test,"- causal"), round(d$mcap_human_readable, 2) , round( d$ecap,2), round((d$m -d$mcap) / d$scap, 2), d$pcap_human_readable,  d$symbol_pcap)
			arr_3 = c( paste(cov_to_test,"- effectual"), round(b$mcap_human_readable, 2) , round( b$ecap,2), round((b$m -b$mcap) / b$scap, 2), b$pcap_human_readable,  b$symbol_pcap)
				
			writeLines("------------------------------------------------------------------------------")
			writeLines(sprintf("Beta-Coefficients for variable no. %s:", i))
			writeLines("\t\t\tEstimate\tStd. Error\tz value\t\tPr(>|z|)")
			writeLines(sprintf("%s\t\t%s\t\t%s\t\t%s\t\t%s %s", arr_1[1], arr_1[2], arr_1[3], arr_1[4], arr_1[5], arr_1[6] ))
			writeLines(sprintf("%s\t\t%s\t\t%s\t\t%s\t\t%s %s", arr_2[1], arr_2[2], arr_2[3], arr_2[4], arr_2[5], arr_2[6] ))
			writeLines(sprintf("%s\t\t%s\t\t%s\t\t%s\t\t%s %s", arr_3[1], arr_3[2], arr_3[3], arr_3[4], arr_3[5], arr_3[6] ))
			
			my_df = as_tibble(rbind(arr_1, arr_2, arr_3))
			colnames(my_df) = arr_names
			my_df$variable = i
			my_df$algorithm = "sr"
			my_df

		}))
	
		writeLines("------------------------------------------------------------------------------")
		writeLines("")
		writeLines("Significance codes: 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
		writeLines(sprintf("Number of Observations: %s", nrow(my_design)))
		writeLines("Link: linear-softmax")
	}
	
	# Plots
	plot(grid.arrange(
		plot_densities(dt$coef_ang_posterior , do_log = F, color="0",  fill = 1:ncol(proportions)) +
		ggtitle("Extrinsic trends - posterior distribution of angular coeff. - Dirichlet regression"),
		plot_densities(bt$coef_ang_posterior , do_log = F, color="0",  fill = 1:ncol(proportions)) +
			ggtitle("Intrinsic trends - posterior distribution of angular coeff. - Beta regression")
	))


	list(fit = fit, coef.mat = coef.mat)

}
