

get_estimates = function(.data, lev,X) {
	
	.data %>% 
		filter(level ==lev) %>%
		filter(.variable %>% is.na %>% `!`) %>%
		select(level, `Cell type category`, draws) %>% 
		mutate(regression = map(draws,
														~ .x %>%
															group_by(A) %>%
															summarise(.median = median(.value), .sd = sd(.value)) %>% 
															#tidybayes::median_qi(.width = credible_interval) %>%
													
															left_join(tibble(A=1:ncol(X), A_name = colnames(X)) ,  by = "A") %>% 
															select(-A) %>% 
														
															pivot_wider(
																names_from = A_name,
																values_from = c(median, sd)
															))) %>% 
		select(-draws) %>% 
		unnest(regression)
	
}
