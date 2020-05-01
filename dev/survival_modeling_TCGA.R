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


# All cancers PFI plot
PFI =  dir(my_dir, full.names = T) %>%
	grep("_Primary_Tumor", ., value = T) %>%
	map_dfr(~.x %>% read_csv %>% distinct(cases_0_project_disease_type, sample))  %>%
	tidyr::extract(sample, into = "sample", regex = "([a-zA-Z0-9]+-[a-zA-Z0-9]+-[a-zA-Z0-9]+)") %>%
	left_join(
		read_csv("dev/survival_TCGA_curated.csv") %>%
			select(bcr_patient_barcode, type, PFI.2, PFI.time.2) %>%
			mutate(PFI.2 = ifelse(PFI.2 == "#N/A", NA, PFI.2)) %>%
			mutate(PFI.time.2 = ifelse(PFI.time.2 == "#N/A", NA, PFI.time.2)) %>%
			mutate(PFI.2 = as.integer(PFI.2), PFI.time.2 = as.integer(PFI.time.2)),
		by = c("sample" = "bcr_patient_barcode")
	)

PFI %>%
	saveRDS("dev/PFI_all_cancers.rds", compress = "gzip")

PFI %>%
	ggplot(aes(PFI.time.2, fill=factor(PFI.2))) + geom_histogram() + facet_wrap(~cases_0_project_disease_type) + scale_x_log10()

PFI_dead = brm(
	PFI.time.2 ~ 0 + cases_0_project_disease_type , 
	data = PFI %>% filter(PFI.2 == 1) %>% filter(PFI.time.2 > 0), 
	family = "gamma", 
	cores = 4, 
	iter = 500, 
	prior = set_prior("normal(mu_d, sigma_d)"), 
	stanvars = stanvar(scode = "real mu_d; real<lower=0> sigma_d;", block = "parameters")
)

PFI_censored = brm(
	PFI.time.2 | cens(alive) ~ 1, # 0 + cases_0_project_disease_type , 
	data = PFI %>% mutate(alive = PFI.2 == 0) %>% filter(PFI.time.2 > 0), # %>% inner_join( (.) %>% distinct(cases_0_project_disease_type) %>% slice(1:5)), 
	family = "gamma", 
	cores = 4
	# , 
	# iter = 500, 
	# prior = 
	# 	set_prior("normal(mu_d, sigma_d)") + 
	# 	set_prior("target += student_t_lpdf(sigma_d | 3, 0, 2)", check = FALSE) +
	# 	set_prior("target += student_t_lpdf(mu_d | 3, 0, 2)", check = FALSE) , 
	# stanvars = stanvar(scode = "real mu_d; real<lower=0> sigma_d;", block = "parameters")
)
