library(tidyverse)
library(ARMET)
# library(furrr)
library(tidybulk)
library(dendextend)
library(RColorBrewer)
library(nanny)
library(tidyHeatmap)
library(ggrepel)
# plan(multicore)
# library(doParallel)
# registerDoParallel(20)

library(furrr)
plan(multisession, workers=15)
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
		theme(legend.position="bottom", aspect.ratio = 1)
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
	dir("dev/armet_TCGA_Ju29_lv1_beta_lv2_dir/", pattern = "^armet_", full.names = T) %>%	grep("regression.rds$", ., value = T) %>%
	map_dfr(~ readRDS(.x) %>%
						# Add relative probability of slope != 0
						mutate(alpha2_prob =
									 	map_dbl(
									 		draws_cens, 
									 		~ .x %>% 
									 			distinct() %>%
									 			#filter(A == 2) %>%
									 			ARMET:::draws_to_prob_non_zero()
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
							.value_alpha1, .value_alpha2,.lower_alpha2, .upper_alpha2, .lower_alpha2_cens, .upper_alpha2_cens,
							.value_alpha2_cens,
							alpha2_prob
						)))
					
	)

# hmap_df %>% saveRDS("dev/hmap_df.rds") 
cancer_imm_therapy = c("LUAD", "LUSC", "SKCM", "metastatic_SKCM", "KIRC", "KICH", "KIRP", "HNSC")
hmap_df_annotated =
	hmap_df %>%
	extract(file, c( "cancer_ID", "dummy"), regex = ".*armet_((metastatic_)?[A-Z]+).*") %>%
	select(-dummy) %>%
	
	# Add annotation
	left_join(read_csv("dev/TCGA_supergroups.csv") %>% select(-cancer)) %>% 
	mutate(group = if_else(cancer_ID == "metastatic_SKCM", "melanomas", group)) %>%
	
	mutate(imm_therapy = if_else(cancer_ID %in% cancer_imm_therapy, T, F)) %>%
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
	mutate(.value_alpha1 = .value_alpha1 %>% scale(center = F)) %>%
	heatmap( 
		cancer_ID, `Cell type category`, .value_alpha1, 
		palette_value = circlize::colorRamp2(c(-2, -1, 0, 1, 2), brewer.pal(5, "RdBu")) ,
		.scale = "none"
	) %>%
	add_tile(group, palette = brewer.pal(8, "Set1")) %>%
	add_point(median_DSS) %>%
	add_tile(imm_therapy, palette = c("white", "grey"))

hmap_association = 
	hmap_df_annotated %>%
	unnest(data) %>%
	mutate(group = if_else(group %>% is.na, "other", group)) %>%
	mutate(alpha2_combined = alpha2_prob) %>%
	group_by(level) %>%
	heatmap( 
		cancer_ID, `Cell type category`, alpha2_combined ,
		palette_value = circlize::colorRamp2(c(-1, -0.5, 0, 0.5, 1), brewer.pal(5, "RdBu")),
		.scale = "none"
	) %>%
	add_tile(group, palette = brewer.pal(8, "Set1")) %>%
	add_point(median_DSS) %>%
	add_tile(imm_therapy, palette = c("white", "grey"))

hmap_composition %>% save_pdf("dev/hmap_composition.pdf", width = 183, height = 110, units = "mm")
hmap_association %>% save_pdf("dev/hmap_association.pdf", width = 183, height = 110, units = "mm")

library(ComplexHeatmap)
pdf("dev/tanglegram.pdf")
tanglegram(
	hmap_composition %>% show %>% draw %>% row_dend, 
	hmap_association %>% show %>% draw %>% row_dend,
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
		my_theme + theme(aspect.ratio = 1)
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

produce_KM_curves = function(file_name){
	
	file_name %>%
		
		readRDS() %>%
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
		unnest(.draws) %>%
		
		# Stratify the data
		nanny::nest_subset(
			data = -c(type, `Cell type category`, .draw), 
			.exclude = dead
		) %>%
		
		# Define the split
		mutate(data = map(
			data,
			~ .x %>%
				mutate(med_prop = median(.value_relative)) %>%
				mutate(high = .value_relative > med_prop)
		)) %>%
		
		# Execute model
		mutate(fit = map(data,
										 ~
										 	survival::survfit(
										 		survival::Surv(PFI.time.2, dead) ~ high,
										 		data = .x
										 	))) %>%
		
		# Calculate p-value
		nest(ct_data = -c(`Cell type category`)) %>%
		mutate(pvalue = map(
			ct_data,
			~ survminer::surv_pvalue(.x$fit, data = .x$data, combine = TRUE)
		)) %>%
		
		mutate(avg_pvalue = map_dbl(pvalue, ~ median(.x$pval))) %>%
		
		# Execute plot
		mutate(plot_df = map(
			ct_data,
			~ survminer::ggsurvplot(
				fit = .x$fit %>% setNames(as.character(1:length(.))) %>% .[1:150],
				.x$data  %>% setNames(as.character(1:length(.))) %>% .[1:150],
				risk.table = FALSE,
				conf.int = F,
				combine = TRUE,
				legend = "none",
				pval = F
			)$plot$data
		)) %>%
		
		mutate(plot = map(
			plot_df,
			~ .x %>%
				separate(strata, c("id", "category"), sep = "::", remove = F) %>%
				arrange(surv, id, category) %>%
				ggplot(aes(time, surv,  color = category)) +
				geom_line(aes(group = strata), alpha = 0.2) +
				my_theme + theme(text = element_text(size = 8), aspect.ratio = 1)
		))
}

produce_KM_curves_abs_values = function(file_name){
	
	file_name %>%
		
		readRDS() %>%
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
		#unnest(.draws) %>%
		
		# Stratify the data
		nanny::nest_subset(
			data = -c(type, `Cell type category`, level), 
			.exclude = dead
		) %>%
		
		# Define the split
		mutate(data = map(
			data,
			~ .x %>%
				mutate(med_prop = median(proportion)) %>%
				mutate(high = proportion > med_prop)
		)) %>%
		
		# Execute model
		mutate(fit = map(data,
										 ~
										 	survival::survfit(
										 		survival::Surv(PFI.time.2, dead) ~ high,
										 		data = .x
										 	))) %>%
		
		# Calculate p-value
		mutate(pvalue = map2(
			fit, data,
			~ survminer::surv_pvalue(.x, data = .y)
		)) %>%
		
		# Execute plot
		mutate(plot = map2(
			fit, data,
			~ survminer::ggsurvplot(
				fit = .x,
				.y,
				risk.table = FALSE,
				conf.int = T,
				#legend = "none",
				pval = T, 
				palette = "Set1", 
				ggtheme = my_theme + theme(axis.title.y=element_blank(),
																	 axis.text.y=element_blank(), aspect.ratio = 1)
			) 
		)) 
}


# Kaplan-Meyer curves
km_curves =
	dir("dev/armet_TCGA_Ju29_lv1_beta_lv2_dir/", pattern = "^armet_", full.names = T) %>%	grep("regression.rds$", ., value = T) %>%
	
	# Prepare data
	enframe(value = "file_name") %>%
	
	mutate(km_data = future_map(file_name, ~ .x %>% produce_KM_curves_abs_values() )) %>%
	unnest(km_data)

km_curves %>% saveRDS("dev/km_curves.rds")


library(patchwork)
library(survminer)

km_grid = 
	km_curves %>%
	tidyr::extract(file_name, "type", ".*armet_([A-Z]+)\\..*") %>%
	filter(type %in% cancer_imm_therapy) %>%
	filter(level<=3) %>%
	# Get top per cancer
	unnest(pvalue) %>%
	group_by(type) %>%
	arrange(pval) %>%
	slice(1) %>%
	ungroup() %>%
	filter(pval <0.05) %>%
	
	mutate(plot = 	
				 	pmap(list(plot, type, `Cell type category`), ~	..1 +
				 			 	#	scale_fill_brewer(palette="Set1") +
				 			 	ggtitle(sprintf("%s %s", ..2, ..3)) 
				 	)) %>%
	pull(plot) %>%
	arrange_ggsurvplots(ncol = 6, nrow = 1)


ggsave(
	"dev/KM_curves_ARMET_TCGA.pdf",
	plot = km_grid,
	useDingbats=FALSE,
	units = c("mm"),
	width = 200 ,
	height = 183*0.61,
	limitsize = FALSE
)

# Rank of immunogenicity
(
	hmap_df_annotated %>% 
		unnest(data) %>%
		#	extract(file_name, "type", ".*armet_([A-Z]+)\\..*") %>%
		filter(`Cell type category` == "immune_cell") %>%
		arrange(.value_alpha2_cens %>% desc) %>%
		mutate(cancer_ID = factor(cancer_ID, levels = unique(.$cancer_ID))) %>%
		ggplot(aes(cancer_ID, .value_alpha2_cens)) +
		geom_hline(yintercept = 0, type="dashed", color="grey") +
		geom_errorbar(aes(ymin=.lower_alpha2_cens, ymax=.upper_alpha2_cens), width=0) +
		geom_point(aes( size=.value_alpha1, fill=imm_therapy), shape=21) + 
		coord_flip(ylim = c(-1, 1.3)) +
		scale_fill_manual(values = c( "grey",  "yellow")) +
		my_theme + theme(aspect.ratio = 1)
) %>%
	ggsave(
		"dev/immunogenicity_TCGA.pdf",
		plot = .,
		useDingbats=FALSE,
		units = c("mm"),
		width = 183 ,
		height = 183,
		limitsize = FALSE
	)


# Plot of best ranked
dd =
	dir("dev", pattern = "^armet_", full.names = T) %>%	grep("regression.rds$", ., value = T) %>%
	# Prepare data
	map_dfr(
		~ readRDS(.x) %>% mutate(file=.x))

dd %>%
	extract(file, "type", ".*armet_([A-Z]+)\\..*") %>%
	inner_join(
		hmap_df_annotated %>% 
			unnest(data) %>% 
			arrange(alpha2_prob %>% abs %>% desc, .value_alpha2_cens %>% abs %>% desc) %>%
			select(type = cancer_ID, `Cell type category`) %>% 
			slice(1:6)
	) %>%
	select( `Cell type category`, proportions) %>% 
	unnest(proportions) %>% 
	ggplot(aes(boot::logit(.value_relative), (PFI.time.2))) + 
	geom_errorbar(aes(xmin = boot::logit(.value_relative.lower), xmax = boot::logit(.value_relative.upper))) + 
	geom_point() + 
	facet_wrap(~ type + `Cell type category`, scale="free_x") +
	scale_y_log10()


# P-value histogram
km_curves <- readRDS("dev/km_curves.rds")
km_cibersort  = readRDS("dev/km_cibersort.rds")

set.seed(321)
(
	km_curves %>% 
		select(pvalue) %>% 
		unnest(pvalue) %>%
		select(pval) %>%
		mutate(algorithm="ARMET") %>%
		bind_rows(
			km_cibersort %>% 
				unnest(pvalue) %>%
				select(pval) %>%
				mutate(algorithm="Cibersort")
		) %>%
		ggplot(aes(pval)) +
		geom_histogram(aes(y=..count../sum(..count..))) +
		facet_wrap(~algorithm) +
		my_theme +
		theme(	axis.text.x = element_text(	angle = 0), aspect.ratio = 1)
) %>%
	ggsave(
		"dev/TCGA_KM_curves_pvalue_hist.pdf.pdf",
		plot = .,
		useDingbats=FALSE,
		units = c("mm"),
		width = 183 ,
		height = 183*0.61,
		limitsize = FALSE
	)

km_curves %>%
	unnest(pvalue) %>%
	arrange(pval) %>%
	slice(1:3) %>%
	pull(plot) %>%
	survminer::arrange_ggsurvplots() +
	plot_layout(	guides = "collect", nrow = 1	)  & 
	theme(legend.position = 'bottom') 

km_cibersort %>%
	unnest(pvalue) %>%
	arrange(pval) %>%
	slice(1:3) %>%
	mutate(plot = map2(
		fit, data,
		~ survminer::ggsurvplot(
			fit = .x,
			.y,
			risk.table = FALSE,
			conf.int = T,
			#legend = "none",
			pval = T
		)
	)) %>%
	pull(plot)

# Example fit Cibersort-ARMET

# km_curves %>% 
# 	unnest(pvalue) %>%
# 	arrange(pval) %>% 
# 	slice(1:3) %>% 
# 	unnest(data) %>% 
# 	mutate(alive = !dead) %>%
# 	select(-dead) %>%
# 	mutate(method="ARMET") %>%
# 	rename(.cell_type = `Cell type category`) %>%
# 	bind_rows(
# 		km_cibersort %>% 
# 			arrange(p.value) %>% 
# 			slice(1:3) %>% 
# 			unnest(data) %>%
# 			mutate(method="Cibersort") %>%
# 			rename(proportion = .proportion) 
# 	) %>%
# 	ggplot(aes(proportion, PFI.time.2, color=alive)) + 
# 	geom_point() + 
# 	geom_smooth(method="lm") +
# 	facet_wrap(method~interaction(type, .cell_type)) +
# 	scale_y_log10() + 
# 	scale_color_manual(values=c( "#e11f28", "grey")) +
# 	scale_x_continuous(trans="probit") + 
# 	my_theme

base_breaks <- function(n = 10){
	function(x) {
		axisTicks(boot::logit(range(x, na.rm = TRUE)), log = FALSE, n = n)
	}
}


library("functional")
library(scales)
logit <-
	trans_new("logit",
						transform = qlogis,
						inverse = plogis,
						breaks = Compose(qlogis, extended_breaks(), plogis),
						format = scales::label_scientific(digits = 2)
	)

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
													 draw_group = function(self, data, ..., draw_quantiles = NULL) {
													 	data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
													 	grp <- data[1, "group"]
													 	newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
													 	newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
													 	newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
													 	
													 	if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
													 		stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
													 																							1))
													 		quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
													 		aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
													 		aesthetics$alpha <- rep(1, nrow(quantiles))
													 		both <- cbind(quantiles, aesthetics)
													 		quantile_grob <- GeomPath$draw_panel(both, ...)
													 		ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
													 	}
													 	else {
													 		ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
													 	}
													 })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
															draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
															show.legend = NA, inherit.aes = TRUE) {
	layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
				position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
				params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}


p1 = (
	km_cibersort %>% 
		arrange(p.value) %>% 
		slice(c(1)) %>% 
		unnest(data) %>% 
		ggplot(aes(.proportion_0_corrected, PFI.time.2, color=alive)) + 
		geom_point() + 
		geom_smooth(method="lm") +
		scale_y_log10() + 
		scale_color_manual(values=c( "#e11f28", "grey")) +
		scale_x_continuous(trans = logit) + 
		my_theme+ theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),  axis.text.x = element_text(
			angle = 40,
			hjust = 1,
			vjust = 1
		), aspect.ratio = 1)
) %>% 
	ggExtra::ggMarginal(type = "histogram")

p2 = (
	km_cibersort %>% 
		arrange(p.value) %>% 
		slice(2) %>% 
		unnest(data) %>% 
		ggplot(aes(.proportion_0_corrected, PFI.time.2, color=alive)) + 
		geom_point() + 
		geom_smooth(method="lm") +
		scale_y_log10() + 
		scale_color_manual(values=c( "#e11f28", "grey")) +
		scale_x_continuous(trans = logit) + 
		my_theme+ theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),  axis.text.x = element_text(
			angle = 40,
			hjust = 1,
			vjust = 1
		), aspect.ratio = 1)
) %>% 
	ggExtra::ggMarginal(type = "histogram")

p3 = (
	km_cibersort %>% 
		arrange(p.value) %>% 
		slice(2) %>% 
		unnest(data) %>% 
		ggplot(aes(.proportion_0_corrected, PFI.time.2, color=alive)) + 
		geom_point() + 
		geom_smooth(method="lm") +
		scale_y_log10() + 
		scale_color_manual(values=c( "#e11f28", "grey")) +
		scale_x_continuous(trans = logit) + 
		my_theme+ theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),  axis.text.x = element_text(
			angle = 40,
			hjust = 1,
			vjust = 1
		), aspect.ratio = 1)
) %>% 
	ggExtra::ggMarginal(type = "histogram")

# Example fit ARMET
p4 = (
	
	km_curves %>% 
		unnest(pvalue) %>%
		arrange(pval) %>% 
		slice(1) %>% 
		unnest(data) %>% 
		mutate(alive = !dead) %>%
		ggplot(aes(proportion, PFI.time.2, color=alive)) + 
		geom_point() + 
		geom_smooth(method="lm") +
		scale_y_log10() + 
		scale_color_manual(values=c( "#e11f28", "grey")) +
		scale_x_continuous(trans = logit) + 
		my_theme+ theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),  axis.text.x = element_text(
			angle = 40,
			hjust = 1,
			vjust = 1
		), aspect.ratio = 1)
) %>% 
	ggExtra::ggMarginal(type = "histogram")

p5 =  (
	
	km_curves %>% 
		unnest(pvalue) %>%
		arrange(pval) %>% 
		slice(2) %>% 
		unnest(data) %>% 
		mutate(alive = !dead) %>%
		ggplot(aes(proportion, PFI.time.2, color=alive)) + 
		geom_point() + 
		geom_smooth(method="lm") +
		scale_y_log10() + 
		scale_color_manual(values=c( "#e11f28", "grey")) +
		scale_x_continuous(trans = logit) + 
		my_theme+ theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),  axis.text.x = element_text(
			angle = 40,
			hjust = 1,
			vjust = 1
		), aspect.ratio = 1)
) %>% 
	ggExtra::ggMarginal(type = "histogram")

p6 = 	(
	
	km_curves %>% 
		unnest(pvalue) %>%
		arrange(pval) %>% 
		slice(3) %>% 
		unnest(data) %>% 
		mutate(alive = !dead) %>%
		ggplot(aes(proportion, PFI.time.2, color=alive)) + 
		geom_point() + 
		geom_smooth(method="lm") +
		scale_y_log10() + 
		scale_color_manual(values=c( "#e11f28", "grey")) +
		scale_x_continuous(trans = logit) + 
		my_theme+ theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),  axis.text.x = element_text(
			angle = 40,
			hjust = 1,
			vjust = 1
		), aspect.ratio = 1)
) %>% 
	ggExtra::ggMarginal(type = "histogram")


plot_grid(p1, p2, p3, p4, p5, p6, nrow = 2 ) %>%
	ggsave(
		"dev/example_top_KM_regression_cibersort_ARMET.pdf",
		plot = .,
		useDingbats=FALSE,
		units = c("mm"),
		width = 183/4*3 ,
		height = 150,
		limitsize = FALSE
	)

# Pots composition vs association
data_4_scatter = 
	dir("dev/armet_TCGA_Ju29_lv1_beta_lv2_dir/", pattern = "^armet_", full.names = T) %>%	grep("regression.rds$", ., value = T) %>%
	map_dfr(~ readRDS(.x) %>%
						# Add relative probability of slope != 0
						mutate(alpha2_prob =
									 	map_dbl(
									 		draws_cens, 
									 		~ .x %>% 
									 			distinct() %>%
									 			#filter(A == 2) %>%
									 			ARMET:::draws_to_prob_non_zero()
									 	)) %>%
						
						mutate(file=.x) %>%
						nest(data = -file)
	) %>%
	unnest(data)  %>% 
	mutate(log_avg_proportion = map_dbl(proportions, ~ .x %>% pull(.value) %>% boot::logit() %>% mean %>% boot::inv.logit() )) %>%
	extract(file, c( "cancer_ID", "dummy"), regex = ".*armet_((metastatic_)?[A-Z]+).*") %>%
	mutate(imm_therapy = if_else(cancer_ID %in% cancer_imm_therapy, T, F)) %>%
	mutate(imm_therapy = if_else(cancer_ID %in% cancer_imm_therapy, T, F)) %>%
	left_join(
		read_csv("dev/survival_TCGA_curated.csv") %>% 
			filter(DSS_cr==1) %>%
			group_by(type) %>%
			summarise(log_mean_exp_DSS = mean(log1p(as.numeric(DSS.time.cr)), na.rm=T) %>% exp) %>% 
			arrange(log_mean_exp_DSS),
		by=c("cancer_ID" = "type")
	)


data_4_scatter %>%
	filter(`Cell type category` == "immune_cell") %>%
	ggplot(aes(log_avg_proportion, alpha2_prob, color=log(log_mean_exp_DSS))) +
	geom_point() +
	scale_x_continuous(trans = logit) +
	scale_color_distiller(palette = "Spectral")


# Best cell type
data_4_violin = 
	data_4_scatter %>%
	select(cancer_ID, level, `Cell type category`, log_avg_proportion, alpha2_prob, imm_therapy, log_mean_exp_DSS) %>%
	nest(data = -c(`Cell type category`, level)) %>%
	mutate(mean_prob = map_dbl(data, ~.x%>%pull(alpha2_prob) %>% mean)) %>%
	mutate(sd_prob = map_dbl(data, ~.x%>%pull(alpha2_prob) %>% sd)) %>%
	nest(lv_data = -level) %>%
	mutate(lv_data = map(lv_data, ~.x %>% mutate(sd_prob = scale(sd_prob)))) %>%
	unnest(lv_data) %>%
	arrange(level, mean_prob %>% desc) %>%
	mutate(`Cell type category` = factor(`Cell type category`, levels = unique(.$`Cell type category`))) %>%
	unnest(data) 

p1 = 
	data_4_violin %>%
	ggplot(aes(`Cell type category`, alpha2_prob)) + 
	geom_hline(yintercept = 0) +
	geom_violin() +
	geom_point(alpha=0.4) +
	stat_summary(fun.y=mean, geom="point", size=2, color="red") +
	facet_grid(.~level, scales = "free_x", space = "free_x")+
	theme_bw() +
	theme(
		panel.border = element_blank(),
		axis.line = element_line(),
		panel.grid.major = element_line(size = 0.2),
		panel.grid.minor = element_line(size = 0.1),
		text = element_text(size = 10),
		legend.position = "bottom",
		axis.text.x = element_blank(),
		axis.ticks.x = element_blank(),
		axis.title.x = element_blank(),
		strip.background = element_blank()
	)


p2 = 
	data_4_violin %>%
	ggplot(aes(`Cell type category`, 1, fill=sd_prob)) +
	geom_tile() +
	facet_grid(.~level, scales = "free_x", space = "free_x")+
	scale_fill_gradient2(low ="#440154FF", mid = "#21908CFF", high="#fefada" )  +
	my_theme + theme(aspect.ratio = 1)
#+
#theme(aspect.ratio = 0.03)

library(patchwork)

p = p1 / p2

ggsave(
	"dev/best_cell_types.pdf",
	plot = p,
	useDingbats=FALSE,
	units = c("mm"),
	width = 120 ,
	height = 183,
	limitsize = FALSE
)

# DE analysis for immunotherapies
type_in_cluster_imm = c("STAD", "KIRP", "DLBC", "KICH", "OV", "ACC", "CESC", "TGCT", "HNSC", "LUSC", "MESO", "BLCA", "mSKCM", "SARC")
type_in_cluster_brain = c("LGG", "GBM", "UVM", "PRAD", "KIRC", "PCPG", "THYM")


(
	data_4_scatter %>%
		dplyr::select(cancer_ID, level, `Cell type category`, log_avg_proportion, alpha2_prob, imm_therapy, log_mean_exp_DSS) %>%
		mutate(cluster = case_when( cancer_ID %in% type_in_cluster_imm ~ 1, cancer_ID %in% type_in_cluster_brain ~ 2, TRUE~3 ) %>% factor ) %>%
		nest(ct_data = -`Cell type category` ) %>%
		mutate(fit = map(ct_data, ~ lm( alpha2_prob ~ 0 + cluster , data=.x))) %>%
		mutate(fit_contrast = map(fit, ~ glht(.x, linfct = matrix(c(1,-1/2,-1/2), nrow = 1)))) %>%
		mutate(pvalue = map_dbl(fit_contrast, ~.x %>% summary() %$% test %$% pvalues %>% as.numeric())) %>%
		mutate(coefficient = map_dbl(fit_contrast, ~.x %>% summary() %$% test %$% coefficients %>% as.numeric())) %>%
		
		# plot
		arrange(pvalue) %>%
		mutate(`Cell type category` = factor(`Cell type category`, levels = .$`Cell type category`)) %>%
		filter(`Cell type category` %in% c("t_CD4", "nk_primed", "macrophage_M2", "t_gamma_delta", "nk_primed_IL2", "t_CD8_naive") %>% `!`) %>%
		slice(1:8) %>%
		unnest(ct_data) %>%
		ggplot(aes(cluster, alpha2_prob, fill=cluster == 1)) + 
		geom_boxplot(outlier.shape = NA) +
		geom_point(size=0.3) +
		facet_wrap(~`Cell type category`, nrow = 2) +
		scale_fill_manual(values = c( "grey",  "yellow")) +
		my_theme + theme(aspect.ratio = 1)
) %>%
	ggsave(
		"dev/immunotherapy_cell_types.pdf",
		plot = .,
		useDingbats=FALSE,
		units = c("mm"),
		width = 63 ,
		height = 183,
		limitsize = FALSE
	)


# Gender analyses
gender = 
	dir("dev/armet_TCGA_Oct2_gender/", pattern = "lv_4_gender_regression.rds", full.names = T) %>%
	map_dfr(~ .x %>% readRDS %>% mutate(file=.x)) 


# Abundance
gender %>%
	
	# Order
	arrange(desc(abs(.value_3))) %>% 
	dplyr::slice(1:10) %>%
	
	# unnest
	select(file, 2:5, prob_non_0_3, proportions) %>%
	unnest(proportions) %>%
	
	# Filter rare
	group_by(type, `Cell type category`) %>%
	mutate(mean_abu = mean(boot::logit(.value))) %>%
	ungroup() %>%
	filter(mean_abu > -10) %>%
	
	# plot
	unite("pair", c(`Cell type category`, type), remove = F) %>%
	mutate(pair = factor(pair, levels = unique(.$pair))) %>%
	
	ggplot(aes(pair, .value, fill=gender)) +
	geom_split_violin(trim = TRUE) +
	geom_boxplot( width = 0.25, notch = FALSE, notchwidth = .4, outlier.shape = NA, coef=0) +
	scale_y_continuous(trans=logit) +
	my_theme + theme(aspect.ratio = 1)



# Gender association
gender_slope_df = 
	gender %>%
	
	# Order
	dplyr::arrange(desc(abs(prob_non_0_4))) %>% 
	
	# unnest
	select(file, 2:5, prob_non_0_2, .value_4, prob_non_0_4, proportions) %>%
	
	mutate(proportions = map(
		proportions,
		~ .x %>%
			mutate(med_prop = median(.value_relative)) %>%
			mutate(high = .value_relative > med_prop) %>%
			mutate(dead = !alive) %>%
			unite(col = "high_gender", c(high, gender))
	)) %>%
	
	# Execute model
	mutate(fit = map(proportions,
									 ~
									 	survival::survfit(
									 		survival::Surv(PFI.time.2, dead) ~ high_gender,
									 		data = .x
									 	))) %>%
	
	# Calculate p-value
	mutate(pvalue = map2(
		proportions, fit,
		~ survminer::surv_pvalue(.y, data = .x)
	)) %>%
	
	# Execute plot
	mutate(plot_df =  map2(
		proportions, fit,
		~ survminer::ggsurvplot(
			fit = .y,
			.x ,
			risk.table = FALSE,
			conf.int = T,
			combine = TRUE,
			legend = "bottom",
			pval = F
		)
	)) 

gender_slope_df %>%
	pull(plot_df) %>%
	.[[1]]



plot_km_type_cell_division = function(.data, type, cell, is_high, my_gender){
	
	my_color = case_when(
		is_high & my_gender == "female" ~ "#AF2020",
		!is_high & my_gender == "female" ~ "#377EB8",
		is_high & my_gender == "male" ~ "#7712A0",
		!is_high & my_gender == "male" ~ "#12700D",
		
	)
	
	.data %>%
		
		filter(grepl(type, file)) %>%
		filter(`Cell type category` == cell) %>%
		# Order
		dplyr::arrange(desc(abs(prob_non_0_4))) %>% 
		
		# unnest
		select(file, 2:5, prob_non_0_2, .value_4, prob_non_0_4, proportions) %>%
		
		mutate(proportions = map(
			proportions,
			~ .x %>%
				mutate(med_prop = median(.value_relative)) %>%
				mutate(high = .value_relative > med_prop) %>%
				mutate(dead = !alive) %>%
				mutate(high_gender = high == is_high & gender == my_gender)
		)) %>%
		
		# Execute model
		mutate(fit = map(proportions,
										 ~
										 	survival::survfit(
										 		survival::Surv(PFI.time.2, dead) ~ high_gender,
										 		data = .x
										 	))) %>%
		
		# Calculate p-value
		mutate(pvalue = map2(
			proportions, fit,
			~ survminer::surv_pvalue(.y, data = .x)
		)) %>%
		
		# Execute plot
		mutate(plot_df =  map2(
			proportions, fit,
			~ survminer::ggsurvplot(
				fit = .y,
				.x ,
				risk.table = FALSE,
				conf.int = T,
				combine = TRUE,
				legend = "bottom",
				pval = F, 
				ggtheme = my_theme + theme(aspect.ratio = 1),
				palette =c("#383838", my_color) ,
				censor.size=2, size = 0.2
			)
		)) %>%
		pull(plot_df) %>%
		.[[1]]
}
plot_km_type_cell_division_thym_immune = function(.data, type, cell){
	.data %>%
		
		filter(grepl(type, file)) %>%
		filter(`Cell type category` == cell) %>%
		# Order
		dplyr::arrange(desc(abs(prob_non_0_4))) %>% 
		
		# unnest
		select(file, 2:5, prob_non_0_2, .value_4, prob_non_0_4, proportions) %>%
		
		mutate(proportions = map(
			proportions,
			~ .x %>%
				mutate(med_prop = median(.value_relative)) %>%
				mutate(high = .value_relative > med_prop) %>%
				mutate(dead = !alive) %>%
				mutate(high_gender = case_when(
					high == FALSE & gender == "male" ~ "FALSE_female",
					high == TRUE & gender == "female" ~ "TRUE_male",
					TRUE ~ "rest"
				))
			
		)) %>%
		
		# Execute model
		mutate(fit = map(proportions,
										 ~
										 	survival::survfit(
										 		survival::Surv(PFI.time.2, dead) ~ high_gender,
										 		data = .x
										 	))) %>%
		
		# Calculate p-value
		mutate(pvalue = map2(
			proportions, fit,
			~ survminer::surv_pvalue(.y, data = .x)
		)) %>%
		
		# Execute plot
		mutate(plot_df =  map2(
			proportions, fit,
			~ survminer::ggsurvplot(
				fit = .y,
				.x ,
				risk.table = FALSE,
				conf.int = T,
				combine = TRUE,
				legend = "bottom",
				pval = F, 
				ggtheme = my_theme + theme(aspect.ratio = 1),
				palette = c( "#377EB8", "#383838", "#7712A0"),
				censor.size=2, size = 0.2
			)
		)) %>%
		pull(plot_df) %>%
		.[[1]]
}

list(
	plot_km_type_cell_division(gender, "KIRC", "mono_derived", TRUE, "female") + ggtitle(paste(c("KIRC", "mono_derived", TRUE, "female"), collapse=" ")),
	plot_km_type_cell_division(gender, "COAD", "epithelial", TRUE, "male") + ggtitle(paste(c("COAD", "epithelial", TRUE, "male"), collapse=" ")),
	plot_km_type_cell_division(gender, "COAD", "fibroblast", FALSE, "male") + ggtitle(paste(c("COAD", "fibroblast", FALSE, "male"), collapse=" ")),
	plot_km_type_cell_division(gender, "KIRC", "fibroblast", TRUE, "male") + ggtitle(paste(c("KIRC", "fibroblast", TRUE, "male"), collapse=" ")),
	plot_km_type_cell_division(gender, "COAD", "immune_cell", FALSE, "female") + ggtitle(paste(c("COAD", "immune_cell", FALSE, "female"), collapse=" ")),
	plot_km_type_cell_division(gender, "HNSC", "immune_cell", TRUE, "male") + ggtitle(paste(c("HNSC", "immune_cell", TRUE, "male"), collapse=" ")),
	plot_km_type_cell_division(gender, "GBM", "mast_cell", TRUE, "female") + ggtitle(paste(c( "GBM", "mast_cell", TRUE, "female"), collapse=" ")),
	plot_km_type_cell_division(gender, "LGG", "mast_cell", FALSE, "female") + ggtitle(paste(c( "LGG", "mast_cell", FALSE, "female"), collapse=" ")),
	plot_km_type_cell_division(gender, "THYM", "epithelial", FALSE, "female") + ggtitle(paste(c("THYM", "epithelial", FALSE, "female"), collapse=" ")),
	plot_km_type_cell_division_thym_immune(gender, "THYM", "immune_cell") + ggtitle(paste(c("THYM", "immune_cell"), collapse=" "))
) %>%
	survminer::arrange_ggsurvplots(ncol=3, nrow=3) %>%
	ggsave(
		"dev/gender_KM.pdf",
		plot = .,
		useDingbats=FALSE,
		units = c("mm"),
		width = 230 ,
		height = 230,
		limitsize = FALSE
	)
library(tidybulk)	

# MA plot
(
	gender %>%
		mutate(type = map_chr(proportions, ~ .x %>% distinct(type) %>% pull(type))) %>%
		# Order
		dplyr::arrange(desc(abs(prob_non_0_4))) %>%
		
		# Setup labels
		mutate(`Cell type category` = if_else(abs(prob_non_0_4) > 0.9, `Cell type category`, "")) %>%
		mutate(type = if_else(abs(prob_non_0_4) > 0.9, type, "")) %>%
		
		ggplot(aes(-.value_alpha4, 1-abs(prob_non_0_4), label=type)) + 
		geom_point(aes(color=`Cell type category`, size = abs(prob_non_0_4)>0.9, alpha=abs(prob_non_0_4)>0.9)) +
		ggrepel::geom_text_repel(size=2.5) +
		#scale_color_brewer(palette="Set1") +
		
		# Custom scales
		scale_y_continuous(trans = "log10_reverse") +
		scale_size_discrete(range = c(1, 2)) +
		scale_color_manual(values=c("black", "#E41A1C" ,"#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")) +
		my_theme + theme(aspect.ratio = 1)
) %>%
	ggsave(
		"dev/gender_MA_plot.pdf",
		plot = .,
		useDingbats=FALSE,
		units = c("mm"),
		width = 89 ,
		height = 100,
		limitsize = FALSE
	)

# PCA f genders
library(plotly)
(
	gender %>%
		mutate(type = map_chr(proportions, ~ .x %>% distinct(type) %>% pull(type))) %>%
		mutate(female = .value_alpha1, male = .value_alpha1+.value_alpha3) %>%
		select(type, `Cell type category`, female, male) %>%
		
		
		gather(sex, value, c(female, male)) %>%
		nanny::reduce_dimensions(c(sex, type), `Cell type category`, value, scale = FALSE, method="PCA", action="get", .dims = 4, ) %>%
		nest(data = -type) %>%
		mutate(dist = map_dbl(data, ~ .x %>% nanny::as_matrix(rownames = sex) %>% dist )) %>%
		unnest(data) %>%
		
		left_join(
			read_csv("dev/TCGA_supergroups.csv") %>% 
				select(-cancer) , 
			by = c("type" = "cancer_ID")
		) %>% 
		mutate(group = if_else(group %>% is.na, "other", group)) %>%
		
		mutate(type_label = if_else(sex == "female", type, "")) %>%
		
		ggplot(aes(PC1, PC2, group=type, shape=sex, label = type_label)) +
		geom_line(alpha=0.8) +
		geom_point(aes(size=dist, color=group), alpha=0.8) + 
		ggrepel::geom_text_repel(size=2, segment.linetype = 2, segment.curvature = -0.1) +
		scale_shape_manual(values =  c(3, 17)) +
		scale_size_continuous(range =  c(1,4)) +
		scale_color_manual(values =  colorRampPalette(brewer.pal(8, "Dark2"))(9)) +
		my_theme + theme(aspect.ratio = 1)
) %>%
	ggsave(
		"dev/gender_PCA_plot.pdf",
		plot = .,
		useDingbats=FALSE,
		units = c("mm"),
		width = 89 ,
		height = 100,
		limitsize = FALSE
	)

# Violin plot
(
	gender %>%
		mutate(type = map_chr(proportions, ~ .x %>% distinct(type) %>% pull(type))) %>%
		mutate(female = .value_alpha1, male = .value_alpha1+.value_alpha3) %>%
		select(type, `Cell type category`, female, male, diff = .value_alpha3) %>%
		
		nest(data = -type) %>%
		mutate(mean_abs_diff = map_dbl(data, ~ .x %>% pull(diff) %>% abs %>% mean )) %>%
		
		
		left_join(
			read_csv("dev/TCGA_supergroups.csv") %>% 
				select(-cancer) , 
			by = c("type" = "cancer_ID")
		) %>% 
		mutate(group = if_else(group %>% is.na, "other", group)) %>%
		
		arrange(mean_abs_diff %>% desc) %>%
		mutate(type = factor(type, levels = unique(.$type))) %>%
		unnest(data) %>%
		
		mutate(cell_type_label = if_else(abs(diff)>3, `Cell type category`, "")) %>%
		
		ggplot(aes(type, diff, label = cell_type_label)) +
		geom_violin() +
		geom_point(aes(color=group), size = 0.7, alpha=0.8) + 
		ggrepel::geom_text_repel(size=2, segment.linetype = 5, segment.curvature = -0.1) +
		scale_color_manual(values =  colorRampPalette(brewer.pal(8, "Dark2"))(9)) +
		my_theme
) %>%
	ggsave(
		"dev/gender_rank_plot.pdf",
		plot = .,
		useDingbats=FALSE,
		units = c("mm"),
		width = 183/100*60 ,
		height = 100,
		limitsize = FALSE
	)

unnest(proportions) %>%
	
	# Filter rare
	group_by(type, `Cell type category`) %>%
	mutate(mean_abu = mean(boot::logit(.value))) %>%
	ungroup() %>%
	filter(mean_abu > -10) %>%
	
	# plot
	unite("pair", c(`Cell type category`, type), remove = F) %>%
	mutate(pair = factor(pair, levels = unique(.$pair))) %>%
	
	ggplot(aes(pair, .value, fill=gender)) +
	geom_split_violin(trim = TRUE) +
	geom_boxplot( width = 0.25, notch = FALSE, notchwidth = .4, outlier.shape = NA, coef=0) +
	scale_y_continuous(trans=logit) +
	my_theme




# Stratify the data
nanny::nest_subset(
	data = -c(type, `Cell type category`, .draw), 
	.exclude = dead
) %>%
	
	# Define the split
	