library(tidyverse)
library(tidyBulk)

# Bhupinder samples
path = "~/third_party_analyses/Bhupinder_oligo_breast"
load(sprintf("%s/results_seurat_together.rda", path))


tt =
	dir(
		sprintf("%s/bulkRNA/alignment_hg38", path),
		pattern = "bam",
		recursive = T,
		full.names = T
	) %>%
		grep("Undetermined", ., invert = T, value = T) %>%
		tidyBulk_SAM_BAM(genome = "hg38")

tt %>%
	filter(transcript %>% is.na %>% `!`) %>%
	scale_abundance() %>%
	ggplot(aes(count_scaled +1, group=sample)) +
	geom_density() +
	scale_x_log10()

