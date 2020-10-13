	library(tidyverse)
	library(magrittr)
	library(purrr)
	library(ARMET)
	
	my_theme =
		theme_bw() +
		theme(
			panel.border = element_blank(),
			axis.line = element_line(),
			panel.grid.major = element_line(size = 0.2),
			panel.grid.minor = element_line(size = 0.1),
			text = element_text(size = 12),
			legend.position = "bottom",
			aspect.ratio = 1,
			strip.background = element_blank(),
			axis.title.x  = element_text(margin = margin(
				t = 10,
				r = 10,
				b = 10,
				l = 10
			)),
			axis.title.y  = element_text(margin = margin(
				t = 10,
				r = 10,
				b = 10,
				l = 10
			)),
			axis.text.x = element_text(angle = 90, hjust = 1)
		)
	
	
	my_dir = "dev/test_simulation"
	
	
	roc_df = 
		
		dir(my_dir, "asses", full.names = TRUE) %>%
		map_dfr(~ readRDS(.x)) %>%
		
		# Calculate ROC
		mutate(abs_slope = abs(slope)) %>%
		nest(data = -c(foreignProp, S, abs_slope, method, CI)) %>%
		mutate(real_negative = map_dbl(data, ~ .x %>% filter(alpha_2==0) %>% nrow)) %>%
		mutate(FP = map_dbl(data, ~ .x %>% filter(fp) %>% nrow)) %>%
		mutate(real_positive = map_dbl(data, ~ .x %>% filter(alpha_2!=0) %>% nrow )) %>%
		mutate(TP = map_dbl(data, ~ .x$tp %>% sum ) ) %>%
		mutate(TP_rate = TP/real_positive, FP_rate = FP/real_negative) %>%
		
		select(-data)
	
	
	
	roc_df %>%
		arrange(FP_rate, TP_rate) %>%
		ggplot(aes(x=FP_rate, y=TP_rate, color=method)) +
		geom_abline(intercept = 0, slope = 1, linetype="dotted", color="grey") +
		geom_line() +
		scale_color_brewer(palette = "Set1") +
		facet_grid(foreignProp ~ abs_slope + S) +
		coord_cartesian(xlim=c(0,0.1), ylim=c(0,1)) +
		my_theme

ggsave(filename = "dev/simulation_benchmark.png")



