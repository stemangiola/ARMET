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
registerDoParallel(20)

library(furrr)
plan(multisession, workers=10)
options(future.globals.maxSize = 50068 * 1024^2)

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

# # Test differential composition
# load("dev/armet_MESO.tcga.harmonized.counts.allgenes.rds.rda")
# dc1 = test_differential_composition(res) 
# load("dev/armet_PAAD.tcga.harmonized.counts.allgenes.rds.rda")
# dc2 = test_differential_composition(res) 
# load("dev/armet_READ.tcga.harmonized.counts.allgenes.rds.rda")
# dc3 = test_differential_composition(res) 
# load("dev/armet_SARC.tcga.harmonized.counts.allgenes.rds.rda")
# dc4 = test_differential_composition(res) 
# load("dev/armet_SKCM.tcga.harmonized.counts.allgenes.rds.rda")
# dc5 = test_differential_composition(res) 
# load("dev/armet_TGCT.tcga.harmonized.counts.allgenes.rds.rda")
# dc6 = test_differential_composition(res) 
# load("dev/armet_THYM.tcga.harmonized.counts.allgenes.rds.rda")
# dc7 = test_differential_composition(res) 
# load("dev/armet_UCS.tcga.harmonized.counts.allgenes.rds.rda" )
# dc8 = test_differential_composition(res) 

dc = 
	foreach(i =	dir("dev", pattern = "^armet_", full.names = T) %>%	grep("rda$", ., value = T) ) %dopar% {
	load(i)
	test_differential_composition(res) %>% saveRDS(sprintf("%s_regression.rds", i))
}

# Polar
pol_p = 
	dir("dev", pattern = "^armet_", full.names = T) %>%	grep("regression.rds$", ., value = T) %>%
	map(~ {
		.x %>%
			readRDS() %>%
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
	#scale_fill_manual(values = cancer_sig %>% distinct(group, color) %>% arrange(group) %>% pull(color)) +
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
	dir("dev", pattern = "^armet_", full.names = T) %>%	grep("regression.rds$", ., value = T) %>%
	map_dfr(~ readRDS(.x) %>%
			# Add relative probability of slope != 0
			mutate(alpha2_prob =
						 	map_dbl(
						 		draws_cens, 
						 		~ .x %>% 
						 			distinct() %>%
						 			#filter(A == 2) %>%
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
			
			mutate(file=.x) %>%
			nest(data = -file) %>%
			# # Polar plot
			# mutate(polar_plot = 
			# 			 	map(data, ~ .x %>%
			# 			 				plot_polar(size_geom_text = 2) +
			# 			 				theme(axis.title.x=element_blank()) +
			# 			 				ggtitle(.x %>% gsub("dev/armet_", "", .) %>% gsub(".rda", "", ., fixed = T)) )) %>%
			
			# Subset for lowering memory
			mutate(data = map(data, ~ .x %>% 			select(
				.variable, community, level ,  
				`Cell type category`, C, 
				.value_alpha1, .value_alpha2, .value_alpha2_cens,
				alpha2_prob
			)))
		
		)

hmap_df %>% 
	saveRDS("hmap_df.rds") 

hmap_df_annotated =
	hmap_df %>%
	extract(file, "cancer_ID", regex = ".*armet_([A-Z]+).*") %>%
	
	# Add annotation
	left_join(read_csv("dev/TCGA_supergroups.csv") %>% select(-cancer)) %>% 
	mutate(imm_therapy = if_else(cancer_ID %in% c("LUAD", "LUSC", "SKCM", "KIRC", "KICH", "KIRP", "HNSC"), T, F)) %>%
	left_join(
		read_csv("dev/survival_TCGA_curated.csv") %>% 
			filter(DSS_cr==1) %>%
			group_by(type) %>%
			summarise(median_DSS = median(as.numeric(DSS.time.cr), na.rm=T)) %>% 
			arrange(median_DSS),
		by=c("cancer_ID" = "type")
	)

hmap_composition = 
	hmap_df_annotated  %>%
	unnest(data) %>%
	mutate(group = if_else(group %>% is.na, "other", group)) %>%
	group_by(level) %>%
	heatmap( 
		cancer_ID, `Cell type category`, .value_alpha1, 
		annotation = c(group, imm_therapy, median_DSS),
		type = c("tile", "tile", "point"),
		#palette_discrete = list(unique(cancer_sig$color)),
		palette_value = circlize::colorRamp2(c(-4, -2, 0, 2, 4)/3*2, brewer.pal(5, "RdBu")) , 
		.scale = "column"
	)

hmap_association = 
	hmap_df_annotated %>%
	unnest(data) %>%
	mutate(group = if_else(group %>% is.na, "other", group)) %>%
	mutate(alpha2_combined = .value_alpha2_cens * abs(alpha2_prob)) %>%
	group_by(level) %>%
	heatmap( 
		cancer_ID, `Cell type category`, alpha2_combined ,
		annotation = c(group, imm_therapy, median_DSS),
		type = c("tile", "tile", "point"),
		#palette_discrete = list(unique(cancer_sig$color)),
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
	#	nest(data = - cancer_ID) %>%
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

print_and_return = function(.data, string){
	print(string)
	.data
}



# Kaplan-Meyer curves
km_curves =
	dir("dev", pattern = "^armet_", full.names = T) %>%	grep("regression.rds$", ., value = T) %>%
	
	# Prepare data
	map_dfr(
		~ readRDS(.x) %>%
			select(`Cell type category`, level, proportions) %>%
			unnest(proportions) %>%
			select(
				type,
				`Cell type category`,
				level,
				sample,
				proportion = .value,
				PFI.time.2,
				dead = PFI.2,
				.draws
			) %>% 
			unnest(.draws)
	) %>%
	print_and_return(".") %>%
	
	# Stratify the data
	nanny::nest_subset(data = -c(type, `Cell type category`, .draw)) %>%
	print_and_return(".") %>%
	
	# Define the split
	mutate(data = map(data, ~ .x %>%
											mutate(med_prop = median(.value_relative)) %>%
											mutate(high = .value_relative > med_prop)
											)) %>%
	# Execute model
	mutate(fit = future_map(data,
												 ~
												 	survival::survfit(
												 	survival::Surv(PFI.time.2, dead) ~ high,
												 	data = .x
												 )
								)) %>%
	print_and_return(".") %>%
	
	# Calculate p-value
	nest(ct_data = -c(type, `Cell type category`)) %>%
	mutate(pvalue = future_map(ct_data, ~ survminer::surv_pvalue(.x$fit, data = .x$data, combine=TRUE	)	)  ) %>%
	print_and_return(".") %>%
	
	mutate(avg_pvalue = map_dbl(pvalue, ~ mean(.x$pval) )) %>%

	# Execute plot
	mutate(plot_df = future_map(
		ct_data,
		~ survminer::ggsurvplot(
			fit = .x$fit %>% setNames(as.character(1:length(.))) %>% .[1:500],
			.x$data  %>% setNames(as.character(1:length(.))) %>% .[1:500],
			risk.table = FALSE,
			conf.int = F,
			combine = TRUE,
			legend = "none",
			pval = F
		)$plot$data
	)) %>%
	print_and_return(".") %>%
	
	mutate(plot = future_map(
		plot_df, ~ .x %>%
			separate(strata, c("id", "category"), sep="::", remove = F) %>% 
			arrange( surv, id, category) %>%
			ggplot(aes(time, surv,  color=category)) + 
			geom_line(aes(group=strata), alpha=0.2) +
			my_theme + theme(text = element_text(size=8))
	))



(
	hmap_df %>% 
	unnest(data) %>%
	extract( col = "file",into =  "type", regex = ".+armet_([A-Z]+)\\..*") %>%
	select(type, `Cell type category`, level, alpha2_prob, .value_alpha2_cens) %>%
	left_join(km_curves, by = c("type", "Cell type category", "level")) %>%
	ggplot(aes(pval, .value_alpha2_cens * abs(alpha2_prob), label =type, ct = `Cell type category`)) +
	geom_point()
) %>% plotly::ggplotly()

 
dd %>% filter(grepl("t_CD4$", `Cell type category`) & grepl("READ", file)) %>%
	select(proportions) %>% unnest(proportions) %>% ggplot(aes(boot::logit(.value_relative), PFI.time.2)) + geom_point() + scale_y_log10()


	survminer::arrange_ggsurvplots(print = FALSE, ncol = 1, nrow=1) %>%
	ggsave(
		"survival_plot_CAPRA_S_publication.pdf", 
		plot = .,
		device = "pdf",
		useDingbats=FALSE,
		units = c("mm"),
		width = 183 ,
		height = 183
	)



