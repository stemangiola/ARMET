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

which_files = dir("dev", pattern = "^armet_", full.names = T) %>%	grep("rda$", ., value = T) %>%
	setdiff(dir("dev", pattern = "^armet_", full.names = T) %>%	grep("regression.rds$", ., value = T) %>% gsub("_regression.rds$", "", .))

dc = 
	foreach(i =	which_files ) %do% {
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

hmap_df %>% 
	saveRDS("dev/hmap_df.rds") 

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
				my_theme + theme(text = element_text(size = 8))
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
				pval = T
			)
		)) 
}


# Kaplan-Meyer curves
km_curves =
	dir("dev", pattern = "^armet_", full.names = T) %>%	grep("regression.rds$", ., value = T) %>%

	# Prepare data
	enframe(value = "file_name") %>%
	
	mutate(km_data = future_map(file_name, ~ .x %>% produce_KM_curves_abs_values() )) %>%
	unnest(km_data)

km_curves %>% saveRDS("dev/km_curves.rds")


library(patchwork)

km_grid = 
	km_curves %>%
	extract(file_name, "type", ".*armet_([A-Z]+)\\..*") %>%
	
	# Get top per cancer
	group_by(type) %>%
	arrange(avg_pvalue) %>%
	slice(1) %>%
	ungroup() %>%
	arrange(avg_pvalue) %>%
	filter(avg_pvalue<0.05) %>%
	mutate(plot = pmap(list(plot, type, `Cell type category`),
										 ~..1 +
										 	scale_color_brewer(palette = "Set1") +
										 	theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
										 	ggtitle(paste(..2, ..3))
										)) %>%
	pull(plot) %>%
	wrap_plots() +
	plot_layout(	guides = "collect", nrow = 1	)  & 
	theme(legend.position = 'bottom') 

ggsave(
	"KM_curves_ARMET_TCGA.pdf",
	plot = km_grid,
	useDingbats=FALSE,
	units = c("mm"),
	width = 183 ,
	height = 183*0.61,
	limitsize = FALSE
)

# Rank of immunogenicity
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
	coord_flip() +
	scale_fill_manual(values = c( "grey",  "yellow")) +
	my_theme


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
		theme(	axis.text.x = element_text(	angle = 0))
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
