# Modeling TCGA survival
library(brms)

load("dev/obj_for_armet_makeflow.RData")

i = "Glioblastoma Multiforme"
i = args[1]
i = gsub("_", " ", i)

my_clin_curated = 
	my_clin %>%
	mutate(DAYS_TO_BIRTH = as.numeric(DAYS_TO_BIRTH)) %>%
	mutate(DFS_MONTHS = as.numeric(DFS_MONTHS)) %>%
	filter(DFS_MONTHS>0) %>%
	mutate(sample = factor(sample)) %>%
	mutate(DAYS_TO_BIRTH = ifelse(is.na(DAYS_TO_BIRTH), mean(DAYS_TO_BIRTH, na.rm = T), DAYS_TO_BIRTH)) %>%
	group_by(bcr_patient_barcode) %>%
	arrange(DFS_MONTHS %>% desc) %>%
	slice(1) %>%
	ungroup() %>%
	filter(DFS_STATUS != "") %>% 
	mutate(DFS_MONTHS_scaled = scale(DFS_MONTHS, center = F)) %>%
	mutate(alive = DFS_STATUS == "DiseaseFree") 


surv_gamma = brm(DFS_MONTHS | cens(alive) ~ 1 , data = my_clin_curated, family = "gamma", cores = 4, iter = 500)

surv_gamma_scaled = brm(DFS_MONTHS_scaled | cens(alive) ~ 1 , data = my_clin_curated, family = "gamma", cores = 4, iter = 500)

surv_gamma_scaled_dead = brm(DFS_MONTHS_scaled ~ 1 , data = my_clin_curated %>% filter(!alive), family = "gamma", cores = 4, iter = 500)

surv_gamma_scaled$fit %>% tidybayes::gather_draws(b_Intercept, shape) %>% summarise(m = mean(.value), s = sd(.value))
