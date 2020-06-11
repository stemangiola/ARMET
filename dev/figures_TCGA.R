library(tidyverse)
library(ARMET)
# library(furrr)
library(tidybulk)
library(dendextend)
library(RColorBrewer)
library(nanny)
library(tidyHeatmap)
# plan(multicore)
library(doParallel)
registerDoParallel(30)

my_theme =
	theme_bw() +
	theme(
		panel.border = element_blank(),
		axis.line = element_line(),
		panel.grid.major = element_line(size = 0.2),
		panel.grid.minor = element_line(size = 0.1),
		text = element_text(size = 10),
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
	facet_wrap(~signature, scales = "free") +
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


hmap_df = 
	foreach(i =
						dir("dev", pattern = "^armet_", full.names = T) %>%
						grep("rda$", ., value = T) ,
					.combine = bind_rows
				) %do% {
		load(i)
		res %>%
			test_differential_composition() %>% 
			
			# Add relative probability of slope != 0
			mutate(alpha2_prob =
						 	map_dbl(
						 		draws, 
						 		~ .x %>% 
						 			distinct() %>%
						 			filter(A == 2) %>%
						 			mutate(higher = .value > 0, lower = .value < 0) %>% 
						 			count(higher) %>%
						 			spread(higher, n) %>%
						 			
						 			# Create column if does not exist
						 			purrr::when(
						 				("TRUE" %in% colnames(.) %>% `!`) ~ mutate(., `TRUE` = 0),
						 				("FALSE" %in% colnames(.) %>% `!`) ~ mutate(., `FALSE` = 0),
						 				~ (.)
						 			) %>%
						 			
						 			# Smaller probability
						 			mutate(prob = min(`FALSE`, `TRUE`)/sum(`FALSE`, `TRUE`)) %>% 
						 			
						 			# Multiply by 2 and invert 
						 			mutate(prob = 1 - (prob * 2)) %>%
						 			
						 			mutate(prob = ifelse(`FALSE`>`TRUE`, -prob, prob)) %>%
						 			pull(prob)
						 )) %>%
			select(
				.variable, community, level ,  
				`Cell type category`, C, 
				.value_alpha1, .value_alpha2,
				alpha2_prob
			) %>%
			mutate(file=i)
	}

hmap_df %>% saveRDS("hmap_df.rds")

hmap_composition = 
	hmap_df %>%
	extract(file, "cancer_ID", regex = ".*armet_([A-Z]+).*") %>%
	left_join(read_csv("dev/TCGA_supergroups.csv") %>% select(-cancer)) %>% 
	mutate(group = if_else(group %>% is.na, "other", group)) %>%
	group_by(level) %>%
	heatmap( 
		cancer_ID, `Cell type category`, .value_alpha1, 
		annotation = group, 
		palette_discrete = list(unique(cancer_sig$color)),
		palette_value = circlize::colorRamp2(c(-4, -2, 0, 2, 4)/3*2, brewer.pal(5, "RdBu")) , 
		.scale = "column"
	)

hmap_association = 
	hmap_df %>%
	extract(file, "cancer_ID", regex = ".*armet_([A-Z]+).*") %>%
	left_join(read_csv("dev/TCGA_supergroups.csv") %>% select(-cancer)) %>%
	mutate(group = if_else(group %>% is.na, "other", group)) %>%
	group_by(level) %>%
	heatmap(
		cancer_ID, `Cell type category`, alpha2_prob ,
		annotation = group,
		palette_discrete = list(unique(cancer_sig$color)),
		palette_value = circlize::colorRamp2(c(-1, -0.5, 0, 0.5, 1), brewer.pal(5, "RdBu")),
		.scale = "none"
	)

hmap_composition %>% save_pdf("dev/hmap_composition.pdf", width = 183, height = 110, units = "mm")
hmap_association %>% save_pdf("dev/hmap_association.pdf", width = 183, height = 110, units = "mm")

pdf("dev/tanglegram.pdf")
tanglegram(
	hmap_composition %>% draw %>% row_dend, 
	hmap_association %>% draw %>% row_dend,
	highlight_branches_lwd = F, 
)
dev.off()
# Heatmap signature

# Association rank

(
	hmap_df %>%
		extract(file, "cancer_ID", regex = ".*armet_([A-Z]+).*") %>%
		nest(data = - cancer_ID) %>%
		mutate(mean_association = map_dbl(data, ~ mean(abs(.x$alpha2_prob)))) %>%
		arrange(desc(mean_association)) %>%
		mutate(cancer_ID = factor(cancer_ID, levels = .$cancer_ID)) %>%
		unnest(data) %>%
		ggplot(aes(cancer_ID, abs(alpha2_prob), group=1)) +
		geom_point() +
		stat_summary(fun.y=mean, geom="line", color="red", size=2) +
		scale_color_brewer(palette = "Set1") +
		xlab("Cancer type") +
		ylab("Association probability") +
		my_theme
) %>%
	ggsave(
		"dev/rank_for_cell_type_association.pdf",
		plot = .,
		useDingbats=FALSE,
		units = c("mm"),
		width = 89 ,
		height = 100,
		limitsize = FALSE
	)



	

