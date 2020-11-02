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
			#aspect.ratio = 1,
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
	
	
	my_dir = "dev/test_simulation_singleCell"
	
	
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

p1 = 
	roc_df %>%
	filter(FP_rate<0.1) %>%
	nest(data = -c(foreignProp, S, abs_slope, method)) %>%
	mutate(auc = map_dbl(data, ~	DescTools::AUC(.x$FP_rate, .x$TP_rate, from = 0, to = 0.1))) %>%
	mutate(method = factor(method, levels = c("ARMET", "cibersort",  "llsr", "epic"))) %>%
	ggplot(aes(abs_slope, auc, color = method)) +
	geom_density(stat = "identity") +
	facet_grid(foreignProp ~  S) + 
	geom_hline(yintercept = 0.005, linetype = "dotted") +
	scale_color_brewer(palette="Set1") +
	my_theme

mix_base = 
	readRDS("~/PhD/deconvolution/ARMET/dev/test_simulation_singleCell_makeflow_pipeline/PBMC_integrated_curated.rds") %>%
	tidyseurat::filter(sample != "SCP591") %>%
	tidyseurat::filter(cell_type_curated != "dendritic_myeloid")

p2 = 
	mix_base %>%
	#sample_n(100000) %>%
	ggplot(aes(UMAP_1, UMAP_2, color = cell_type_curated )) +
	geom_point(shape=".") +
	scale_color_brewer(palette = "Set1") +
	my_theme

friendly_cols <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC")

mix_base %>%
	sample_n(10000) %>%
	plot_ly(
		x = ~`UMAP_1`,
		y = ~`UMAP_2`,
		z = ~`UMAP_3`,
		color = ~cell_type_curated,
			colors = RColorBrewer ::brewer.pal(5, "Set1") ,
		size = 0.2,
		opacity=0.7
	)


ggsave(
	"dev/test_simulation_singleCell_lv4_big.pdf",
	useDingbats=FALSE,
	units = c("mm"),
	width = 183 ,
	height = 183,
	limitsize = FALSE
)



