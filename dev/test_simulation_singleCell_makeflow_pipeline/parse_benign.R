library(purrr)
library(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(tidyseurat)
library(Seurat)
library(tidysc)
library(ggplot2)

my_dir = "~/PostDoc/oligo_breast/expanded_analyses_with_control"
sample_info = read_csv(sprintf("%s/OMBC - Sheet1.csv", my_dir))

my_theme =
	theme_bw() +
	theme(
		panel.border = element_blank(),
		axis.line = element_line(),
		panel.grid.major = element_line(size = 0.2),
		panel.grid.minor = element_line(size = 0.1),
		text = element_text(size=8),
		legend.position="bottom",
		
		strip.background = element_blank()
	)

# Load Saskia dataset = 1 SAMPLE
load(sprintf("%s/../data/benign_PBMC/saskia_GSE115189/Sce_Dataset1.RData", my_dir))

# add gene symbol
sce_count = sce@assays$data@listData$counts
rownames(sce_count) = tibble(ens = rownames(sce_count)) %>% tidybulk::ensembl_to_symbol(ens) %>% pull(transcript)
sce_count = sce_count[!is.na(rownames(sce_count)), ]
sce_count = sce_count[!duplicated(rownames(sce_count)), ]

saskia_PBMC =
	sce_count %>% 
	CreateSeuratObject() %>% list() %>% 
	tidysc:::create_tt_from_seurat(species = "Human") %>%
	attr("seurat") %>%
	.[[1]] %>%
	tidyseurat::tidy() %>%
	tidyseurat::mutate(sample = "GSE115189") %>%
	tidyseurat::mutate(source_data = "GSE115189")

# Broad SCP345 - 2 SAMPLES
broad_SCP345 =
	read.delim(sprintf("%s/../data/benign_PBMC/broad_SCP345/pbmc_site_bc.scp.expr.txt", my_dir), row.names = 1) %>% 
	as.matrix %>% 
	Seurat::CreateSeuratObject() %>%
	list() %>%
	tidysc:::create_tt_from_seurat(species = "Human") %>%
	attr("seurat") %>%
	.[[1]] %>%
	tidyseurat::tidy() %>%
	#tidyseurat::select(-sample) %>%
	tidyseurat::left_join(
		read_delim(sprintf("%s/../data/benign_PBMC/broad_SCP345/pbmc_cca_final_metadata_lineage1lineage2.txt", my_dir), "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
			setNames(c("cell",	"Lineages", "kmeans",	"Site" , "sample")) %>%
			mutate(cell = gsub("-", ".", cell))
	) %>%
	tidyseurat::filter(sample %>% is.na %>% `!`) %>%
	tidyseurat::mutate(sample = sprintf("SCP345_%s", sample)) %>%
	tidyseurat::mutate(source_data = "SCP345") %>%
	tidyseurat::rename(old_cell_type = Lineages)

# Broad SCP424 - 2 SAMPLES
counts = Matrix::readMM("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/oligo_breast/data/benign_PBMC/broad_SCP424/matrix.mtx.gz") %>% as.matrix
rownames(counts) = read.csv("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/oligo_breast/data/benign_PBMC/broad_SCP424/features.tsv.gz", sep="", header = F) %>% as_tibble() %>% tidyr::extract(V1, into = "symbol", ".+_(.+)") %>% pull(symbol)
colnames(counts) = read.csv("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/oligo_breast/data/benign_PBMC/broad_SCP424/barcodes.tsv.gz", sep="", header = F)[,1]
counts = counts[!duplicated(rownames(counts)),]

broad_SCP424 =
	counts %>%
	Seurat::CreateSeuratObject() %>%
	list() %>%
	tidysc:::create_tt_from_seurat(species = "Human") %>%
	attr("seurat") %>%
	.[[1]] %>%
	tidyseurat::tidy() %>%
	tidyseurat::left_join(
		read_delim(sprintf("%s/../data/benign_PBMC/broad_SCP424/meta_data_broad_PBMC.txt", my_dir), "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
			setNames(c("cell" ,  "nGene" , "nUMI" , "percent.mito" , "Cluster", "CellType" , "Experiment", "Method" )) %>%
			filter(Method %in% c("10x Chromium (v3)", "10x Chromium (v2)"))
	) %>%
	tidyseurat::filter(Method %>% is.na %>% `!`) %>%
	tidyseurat::mutate(sample = sprintf("SCP424_%s", orig.ident)) %>%
	tidyseurat::mutate(source_data = "SCP424") %>%
	tidyseurat::rename(old_cell_type = CellType)

# Broad SCP589 - ALL depressed not healthy
broad_SCP589 =
	read.delim(sprintf("%s/../data/benign_PBMC/broad_SCP589/Human_PBMC_IFNb_Processed_Log_UMI.txt", my_dir), row.names = 1) %>%
	as.matrix() %>%
	Seurat::CreateSeuratObject() %>%
	list() %>%
	tidysc:::create_tt_from_seurat(species = "Human") %>%
	attr("seurat") %>%
	.[[1]] %>%
	tidyseurat::tidy() %>%
	tidyseurat::left_join(
		read_delim(sprintf("%s/../data/benign_PBMC/broad_SCP589/Human_PBMC_IFNb_Metadata.txt", my_dir), "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
			setNames(c("cell" , 	"sample_condition",	"Cell_Types"	,"Condition"	,"Sex"	,"Age"	,"HCV"	,"CanDep"	,"StimDep"	,"SedDep"	,"CocDep",	"AlcDep",	"NicDep",	"Major_Depressive",	"PTSD" )) %>%
			filter(Condition=="Control")
	) %>%
	tidyseurat::rename(sample = sample_condition)

# Broad SCP591 - 1 SAMPLE
broad_SCP591 = 
	read.delim(sprintf("%s/../data/benign_PBMC/broad_SCP591/Human_PBMC_Morphine_Processed_Log_UMI.txt", my_dir), row.names = 1) %>% 
	as.matrix() %>%
	Seurat::CreateSeuratObject() %>%
	list() %>%
	tidysc:::create_tt_from_seurat(species = "Human") %>%
	attr("seurat") %>%
	.[[1]] %>%
	tidyseurat::tidy() %>%
	tidyseurat::left_join(
		read_delim(sprintf("%s/../data/benign_PBMC/broad_SCP591/Human_PBMC_Morphine_Metadata.txt", my_dir), "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
			setNames(c("cell" , 	"Cell_type",	"group" )) %>%
			filter(group=="Untreated") 
	) %>%
	tidyseurat::filter(group %>% is.na %>% `!`) %>%
	tidyseurat::mutate(sample = "SCP591") %>%
	tidyseurat::mutate(source_data = "SCP591") %>%
	tidyseurat::rename(old_cell_type = Cell_type)

# Cell ranger
cellranger_PBMC = 
	c(
		sprintf("%s/../data/benign_PBMC/10x_6K_PBMC_CallRanger_1_1_0/filtered_gene_bc_matrices/filtered_matrices_mex/hg19/", my_dir),
		sprintf("%s/../data/benign_PBMC/10x_8K_PBMC_CallRanger_2_1_0/filtered_gene_bc_matrices/GRCh38/", my_dir),
		sprintf("%s/../data/benign_PBMC/SRR6260181/outs/filtered_feature_bc_matrix/", my_dir),
		sprintf("%s/../data/benign_PBMC/SRR6260182/outs/filtered_feature_bc_matrix/", my_dir),
		sprintf("%s/../data/benign_PBMC/SRR6260183/outs/filtered_feature_bc_matrix/", my_dir),
		sprintf("%s/../data/benign_PBMC/SRR11038995/outs/filtered_feature_bc_matrix/", my_dir),
		sprintf("%s/../data/benign_PBMC/SRR11038989/outs/filtered_feature_bc_matrix/", my_dir),
		sprintf("%s/../data/benign_PBMC/SRR7722939/outs/filtered_feature_bc_matrix/", my_dir),
		sprintf("%s/../data/benign_PBMC/SRR7244582/outs/filtered_feature_bc_matrix/", my_dir)
	) %>%
	tidysc::tidysc_cell_ranger(	species = "Human"	) %>%
	attr("seurat") %>%
	.[[1]] %>%
	tidyseurat::tidy() %>%
	tidyseurat::extract(orig.ident, "sample", "../data/benign_PBMC/([a-zA-Z0-9_]+)/.+") %>%
	tidyseurat::mutate(sample = dplyr::case_when(
		sample == "10x_6K_PBMC_CallRanger_1_1_0" ~ "10x_6K",
		sample == "10x_8K_PBMC_CallRanger_2_1_0" ~ "10x_8K",
		TRUE ~ sample
	))
	

# Merge
saskia_PBMC %>% 
tidyseurat::bind_rows(cellranger_PBMC) %>%
tidyseurat::bind_rows(broad_SCP345) %>%
tidyseurat::bind_rows(broad_SCP424) %>%
tidyseurat::bind_rows(broad_SCP591) %>%
tidyseurat::bind_rows(broad_SCP589) %>%
tidyseurat::filter(sample %>% is.na %>% `!`) %>%
	
	# Process
	tidyseurat::nest(data = -sample) %>%
	tidyseurat::mutate(data = map(
		data,
		~ .x %>% 
			tidysc::scale_abundance() %>%
			tidysc::reduce_dimensions(method="PCA") %>%
			tidysc::cluster_elements() %>%
			tidyseurat::mutate(cluster = seurat_clusters) %>%
			tidysc::deconvolve_cellularity(species="Human") %>%
			tidysc::reduce_dimensions("UMAP") %>%
			tidysc::reduce_dimensions("tSNE") 
	)) %>%
	saveRDS("dev/test_simulation_singleCell/PBMC_nested.rds", compress = "gzip")


# # PCA
# (
# 	PBMC %>%
# 		
# 	# Annotate
# 	tidyseurat::left_join(
# 		sample_info %>% distinct(sample, race), 
# 		by = "sample"
# 	) %>%
# 		
# 	aggregate_cells() %>%
# 	tidybulk::identify_abundant(sample, transcript, abundance_RNA) %>%
# 	tidybulk::scale_abundance(sample, transcript, abundance_RNA) %>%
# 	tidybulk::reduce_dimensions(sample, transcript, abundance_RNA_scaled, method = "PCA", action="get") %>%
# 	ggplot(aes(PC1, PC2, color=race, label=sample)) + 
# 	geom_text() + 
# 	scale_color_brewer(palette="Set1") +
# 	my_theme + 
# 	theme(aspect.ratio=1)
# ) %>%
# 	plotly::ggplotly()


readRDS("dev/test_simulation_singleCell/PBMC_nested.rds") %>% 
	tidyseurat::unnest(data) %>%
	
	# Harmonize the cell types
	tidyseurat::mutate(
		cell_type_curated = case_when(
			
			old_cell_type == "1. T"    & label_blueprint == "B-cells" ~ "b_cell",  
			old_cell_type == "1. T"    & label_blueprint == "CD4+ T-cells" ~ "t_CD4", 
			old_cell_type == "1. T"    & label_blueprint == "CD8+ T-cells" ~ "t_CD8", 
			old_cell_type == "1. T"    & label_blueprint == "Monocytes" ~ "",  
			old_cell_type == "1. T"    & label_blueprint == "NK cells" ~ "natural_killer", 
			old_cell_type == "2. CD14+ Monocyte"   & label_blueprint == "Monocytes" ~ "monocyte", 
			old_cell_type == "3. NK"    & label_blueprint == "B-cells" ~ "", 
			old_cell_type == "3. NK"    & label_blueprint == "CD4+ T-cells" ~ "", 
			old_cell_type == "3. NK"    & label_blueprint == "CD8+ T-cells" ~ "", 
			old_cell_type == "3. NK"    & label_blueprint == "Monocytes" ~ "", 
			old_cell_type == "3. NK"    & label_blueprint == "NK cells" ~ "natural_killer", 
			old_cell_type == "4. Memory B cell"   & label_blueprint == "B-cells" ~ "b_memory", 
			old_cell_type == "4. Memory B cell"   & label_blueprint == "CD4+ T-cells" ~ "", 
			old_cell_type == "4. Memory B cell"   & label_blueprint == "Monocytes" ~ "", 
			old_cell_type == "4. Memory B cell"   & label_blueprint == "NK cells" ~ "", 
			old_cell_type == "5. DC"    & label_blueprint == "Monocytes" ~ "dendritic_myeloid", 
			old_cell_type == "6. CD16+ Monocyte"   & label_blueprint == "Monocytes" ~ "monocyte", 
			old_cell_type == "7. Naive B cell"   & label_blueprint == "B-cells" ~ "b_naive", 
			old_cell_type == "8. Ambiguous / Potential Doublets" & label_blueprint == "B-cells" ~ "", 
			old_cell_type == "8. Ambiguous / Potential Doublets" & label_blueprint == "CD4+ T-cells" ~ "", 
			old_cell_type == "8. Ambiguous / Potential Doublets" & label_blueprint == "CD8+ T-cells" ~ "", 
			old_cell_type == "8. Ambiguous / Potential Doublets" & label_blueprint == "Monocytes" ~ "", 
			old_cell_type == "8. Ambiguous / Potential Doublets" & label_blueprint == "NK cells" ~ "" , 
			old_cell_type == "B cell"    & label_blueprint == "B-cells" ~ "b_cell" , 
			old_cell_type == "B cell"    & label_blueprint == "CD8+ T-cells" ~ "" , 
			old_cell_type == "B cell"    & label_blueprint == "Monocytes" ~ "" , 
			old_cell_type == "B Cells"    & label_blueprint == "B-cells" ~ "" , 
			old_cell_type == "CD14+ monocyte"   & label_blueprint == "CD8+ T-cells" ~ "" , 
			old_cell_type == "CD14+ monocyte"   & label_blueprint == "Monocytes" ~ "monocyte" , 
			old_cell_type == "CD16+ monocyte"   & label_blueprint == "Monocytes" ~ "monocyte" , 
			old_cell_type == "CD4 T Cells"   & label_blueprint == "B-cells"  ~ "" , 
			old_cell_type == "CD4 T Cells"   & label_blueprint == "CD4+ T-cells" ~ "t_CD4" ,
			old_cell_type == "CD4 T Cells"   & label_blueprint == "CD8+ T-cells" ~ "" , 
			old_cell_type == "CD4+ T cell"   & label_blueprint == "B-cells" ~ "" , 
			old_cell_type == "CD4+ T cell"   & label_blueprint == "CD4+ T-cells" ~ "t_CD4" ,
			old_cell_type == "CD4+ T cell"   & label_blueprint == "CD8+ T-cells" ~ "" , 
			old_cell_type == "CD4+ T cell"   & label_blueprint == "Monocytes" ~ "" , 
			old_cell_type == "CD4+ T cell"   & label_blueprint == "NK cells" ~ "" , 
			old_cell_type == "CD8 T Cells"   & label_blueprint == "CD4+ T-cells" ~ "" ,
			old_cell_type == "CD8 T Cells"   & label_blueprint == "CD8+ T-cells" ~ "t_CD8" , 
			old_cell_type == "Cytotoxic T cell"   & label_blueprint == "CD4+ T-cells" ~ "" , 
			old_cell_type == "Cytotoxic T cell"   & label_blueprint == "CD8+ T-cells" ~ "t_CD8" ,
			old_cell_type == "Cytotoxic T cell"   & label_blueprint == "NK cells" ~ "" ,
			old_cell_type == "Dendritic cell"   & label_blueprint == "CD4+ T-cells" ~ "" , 
			old_cell_type == "Dendritic cell"   & label_blueprint == "Monocytes" ~ "dendritic_myeloid" ,
			old_cell_type == "Megakaryocyte"   & label_blueprint == "CD8+ T-cells" ~ "" ,
			old_cell_type == "Monocytes"    & label_blueprint == "Macrophages" ~ "" ,
			old_cell_type == "Natural killer cell"  & label_blueprint == "CD8+ T-cells" ~ "" , 
			old_cell_type == "Natural killer cell"  & label_blueprint == "NK cells" ~ "natural_killer" , 
			old_cell_type == "NK Cells"    & label_blueprint == "CD8+ T-cells" ~ "" , 
			old_cell_type == "Plasmacytoid dendritic cell" & label_blueprint == "B-cells" ~ "" , 
			old_cell_type == "Plasmacytoid dendritic cell" & label_blueprint == "Monocytes" ~ "" , 
			old_cell_type %>% is.na    & label_blueprint == "Adipocytes" ~ "adipocyte" ,
			old_cell_type %>% is.na    & label_blueprint == "B-cells" ~ "b_cell" ,
			old_cell_type %>% is.na    & label_blueprint == "CD4+ T-cells" ~ "t_CD4" ,
			old_cell_type %>% is.na    & label_blueprint == "CD8+ T-cells" ~ "t_CD8" ,
			old_cell_type %>% is.na    & label_blueprint == "Monocytes" ~ "monocyte" ,
			old_cell_type %>% is.na    & label_blueprint == "NK cells" ~ "natural_killer" ,
			TRUE ~ "THIS SHOULD NOT EXIST"
		)
	) %>% 
	filter(cell_type_curated != "") %>%
	saveRDS("dev/test_simulation_singleCell/PBMC_curated.rds", compress = "gzip")

readRDS("dev/test_simulation_singleCell/PBMC_curated.rds") %>%
	
	tidyseurat::nest(data = -cell_type_curated) %>%
	filter(cell_type_curated != "") %>%
	tidyseurat::mutate(data = map(data, ~ .x %>% tidysc::aggregate_cells(sample))) %>%
	unnest(data) %>%
	
	# Prepare the scaled bulk
	tidyr::unite("sample_ct", c(sample, cell_type_curated), remove = F) %>%
	select(sample_ct, cell_type_curated, transcript, abundance_RNA, sample) %>%
	tidybulk::identify_abundant(sample_ct, transcript, abundance_RNA) %>%
	tidybulk::scale_abundance(sample_ct, transcript, abundance_RNA) %>%
	
	# Save
	saveRDS("dev/test_simulation_singleCell/PBMC_sudo_bulk.rds", compress = "gzip")



readRDS("dev/test_simulation_singleCell/PBMC_sudo_bulk.rds") %>%
	tidybulk::reduce_dimensions(sample_ct, transcript, abundance_RNA_scaled, method="PCA", action="get", scale=F) %>%
	ggplot(aes(PC1, PC2, color=cell_type_curated)) +
	geom_point()
	
