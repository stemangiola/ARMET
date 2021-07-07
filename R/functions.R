



#' test_differential_composition
#' 
#' @description This function performs statistical testing on the fit object
#' 
#' @param .data A tibble
#' @param credible_interval A double
#' @param cluster_CI A double
#' 
#' @export
get_estimates = 
	function(.data, level, credible_interval = 0.90, cluster_CI = 0.55, relative = TRUE) {
		
	
		.d = 
			.data %>%
			filter(.variable %>% is.na %>% `!`) %>%
			extract_CI(credible_interval = credible_interval)
	
		dx = 
			.d %>%
			
			filter(level ==!!level) %>%
			mutate(fold_change_ancestor = 0) %>%
			mutate(significant = ((.lower_alpha2  - 0) * (.upper_alpha2  - 0)) > 0) %>%
			mutate(fold_change  = ifelse(significant, .value_alpha2 , 0))
		
		dx %>% 
			select(
				level, `Cell type category`, .value_alpha1, .value_alpha2, .lower_alpha1,
				.lower_alpha2, .upper_alpha1, .upper_alpha2, significant
			)
		
	}



