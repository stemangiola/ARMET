library(tidyverse)
library(tidyBulk)

# Bhupinder samples
path = "~/PhD/deconvolution/angela_CyTOF_validation/"

tt =
	dir(
		sprintf("%s", path),
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

