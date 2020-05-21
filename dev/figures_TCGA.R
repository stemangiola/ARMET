library(tidyverse)
library(ARMET)
library(furrr)
library(tidybulk)
library(dendextend)
library(RColorBrewer)
library(nanny)
plan(multicore)

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
		axis.text.x = element_text(
			angle = 90,
			hjust = 1,
			vjust = 0.5
		),
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
		))
	)

# Polar
pol_p = 
	dir("dev", pattern = "^armet_", full.names = T) %>%
	grep(".rda", ., fixed = T, value = T) %>%
	future_map(~ {
		load(.x)
		res %>%
			test_differential_composition %>%
			plot_polar(size_geom_text = 2) +
			theme(axis.title.x=element_blank()) +
			ggtitle(.x %>% gsub("dev/armet_", "", .) %>% gsub(".rda", "", ., fixed = T)) 
	}) 

pol_p %>%
	cowplot::plot_grid(plotlist = ., align = "v",  axis="b", rel_widths = 1 ) %>%
	ggsave(
		"dev/landscape_TCGA_polar.pdf",
		plot = .,
		useDingbats=FALSE,
		units = c("mm"),
		width = 347 ,
		height = 347,
		limitsize = FALSE
	)


# PCA
cancer_sig =
	dir("dev", pattern = "^armet_", full.names = T) %>%
	#.[1:3] %>%
	grep("rda", ., fixed = T, value = T) %>%
	future_map_dfr(~ {
		load(.x)
		res %>%
			get_signatures %>%
			mutate(cancer = .x %>% gsub("dev/armet_", "", .) %>% gsub(".rda", "", ., fixed = T))
	}) %>%
	extract(col = cancer, into = "cancer_ID", regex = "([A-Z]+)\\.", remove = F) %>%
	
	# Add info cancer group
	left_join(read_csv("dev/TCGA_supergroups.csv") %>% select(-cancer)) %>%
	
	# Give color
	nest(data = -group) %>%
	arrange(group) %>%
	mutate(color = colorRampPalette(brewer.pal(9, "Set1"))(n())) %>%
	unnest(data) %>%
	
	distinct()
	
saveRDS(cancer_sig, "dev/cancer_sig.rds", compress="gzip")

(
	bind_rows(
	cancer_sig %>%
		filter(A ==1)	%>%
		unite(col="ct_pair", c(`Cell type category_1`, `Cell type category_2`)) %>%
		#	nanny::reduce_dimensions(cancer_ID, c(`Cell type category_1`, `Cell type category_2`), prob, method = "PCA")
		nanny::reduce_dimensions(cancer_ID, ct_pair, prob, method = "PCA") %>%
		nanny::subset(cancer_ID) %>%
		mutate(signature = "Tissue composition"),
	cancer_sig %>%
		filter(A ==2)	%>%
		unite(col="ct_pair", c(`Cell type category_1`, `Cell type category_2`)) %>%
		#	nanny::reduce_dimensions(cancer_ID, c(`Cell type category_1`, `Cell type category_2`), prob, method = "PCA")
		nanny::reduce_dimensions(cancer_ID, ct_pair, prob, method = "PCA") %>%
		nanny::subset(cancer_ID) %>%
		mutate(signature = "Cell-type/Survival")
	) %>%
	mutate(signature = factor(signature, levels=c("Tissue composition", "Cell-type/Survival"))) %>%
	ggplot(aes(PC1, PC2, label=cancer_ID)) +
	geom_point(aes(fill = group), shape = 21, size=2) +
	scale_fill_manual(values = cancer_sig %>% distinct(group, color) %>% arrange(group) %>% pull(color)) +
	ggrepel::geom_text_repel(segment.alpha = 0.5, size=1) +
	facet_wrap(~signature) +
	my_theme +
	theme(legend.position="bottom")
) %>%
	ggsave(
		"dev/PCA_TCGA_groups.pdf",
		plot = .,
		useDingbats=FALSE,
		units = c("mm"),
		width = 183 ,
		limitsize = FALSE
	)

# Heatmap 1
library(doParallel)
registerDoParallel(30)

x = foreach(i = dir("dev/armet_TCGA_may19_before_log_days/", pattern = "^armet_", full.names = T), .combine = bind_rows) %dopar% {
	load(i)
	res %>% test_differential_composition() %>% select(.variable, community, level ,   `Cell type category`, C, .value_alpha1) %>%
		mutate(file=i)
}

x %>%
	extract(file, "cancer_ID", regex = ".*armet_([A-Z]+).*") %>%
	left_join(read_csv("dev/TCGA_supergroups.csv") %>% select(-cancer)) %>% 
	group_by(level) %>%
	heatmap( `Cell type category`, cancer_ID, .value_alpha1, annotation = group, palette_discrete = list(unique(cancer_sig$color)))

# Heatmap 2
library(tidyHeatmap)
cancer_sig %>%
	filter(A == 1) %>%
	nest(data = -c(cancer_ID, cancer, color)) %>%
	mutate(data = map(data, ~ .x %>%
											lower_triangular(`Cell type category_1`, `Cell type category_2`, prob)
											)) %>%
	unnest(data) %>%
	unite("ct", c(`Cell type category_1`, `Cell type category_2`)) %>%
	group_by(level) %>%
	heatmap( ct, cancer_ID, prob, annotation = color, palette_discrete = list(unique(.$color)))



# Dendrograms all features
pdf("dev/TCGA_dendrogram.pdf", useDingbats = F)
cancer_sig %>%
	
	nest(data = -A) %>%
	mutate(
		dendro =
			map(
				data,
				~ .x %>%
					select(-level) %>%
					pivot_wider(
						names_from = c(.variable, `Cell type category_1`, `Cell type category_2`),
						values_from = prob, 
					) %>%
					as_matrix(rownames = "cancer_ID") %>%
					dist() %>%
					hclust() %>%
					as.dendrogram()
			)
	) %>%
	mutate(dendro = map2(
		dendro, data,
		~.x %>% 
			set("nodes_cex", 2) %>%
			set("nodes_pch", 19)  %>% 
			set("leaves_col",  .y %>% distinct(cancer_ID, color) %>% arrange(match(cancer_ID, .x %>% labels))  %>% pull(color)  ) )) %>%
	pull(dendro) %>%
	
	{ tanglegram((.)[[1]], (.)[[2]], highlight_branches_lwd = F) }
dev.off()

# PCA
# pdf("dev/TCGA_dendrogram.pdf", useDingbats = F)
cancer_sig %>%
	nest(data = -A) %>%
	mutate(
		dendro =
			map(
				data,
				~ .x %>%
					unite(col="ct_pair", c(`Cell type category_1`, `Cell type category_2`)) %>%
					nanny::reduce_dimensions(ct_pair, cancer_ID, prob, method = "PCA", scale=F) %>%
					attr("internals") %$%
					PCA %>%
					tcR::pca2euclid(.num.comps = 5) %>%
					as.dist %>%
					hclust() %>%
					as.dendrogram()
			)
	)%>%
	mutate(dendro = map2(
		dendro, data,
		~.x %>% 
			set("nodes_cex", 2) %>%
			set("nodes_pch", 19)  %>% 
			set("leaves_col",  .y %>% distinct(cancer_ID, color) %>% arrange(match(cancer_ID, .x %>% labels))  %>% pull(color)  ) )) %>%
	pull(dendro) %>%
	
	{ tanglegram((.)[[1]], (.)[[2]], highlight_branches_lwd = F) } 
# dev.off()

# Legend
(cancer_sig %>%
		ggplot(aes(prob, level, color = group)) + geom_point() + scale_color_manual(values = cancer_sig %>% distinct(group, color) %>% arrange(group) %>% pull(color))
)%>%
	ggsave(
		"dev/legend_groups.pdf",
		plot = .,
		useDingbats=FALSE,
		units = c("mm"),
		width = 183 ,
		height = 183,
		limitsize = FALSE
	)
