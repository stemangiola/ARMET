library(tidyverse)
library(tidybulk)
library(furrr)
plan(multisession, workers=10)

# Cibersort run

# res_cibersort = 
#   dir("../TCGA_harmonised/", "\\.rds$") %>%
#   enframe(name = NULL, value = "file") %>%
# 
#   mutate(
#     inference =
#       future_map(file, ~ .x %>%
#           prepare_TCGA_input("../") %>%
#           mutate(dead = !alive) %>%
#           tidybulk::test_differential_cellularity(survival::Surv(PFI.time.2, dead) ~ .,
#                                                   sample,
#                                                   transcript,
#                                                   count)
#       )
#   )
# 
# 
#  
# res_cibersort %>% saveRDS("dev/res_cibersort.rds", compress = "gzip")

res_cibersort = readRDS("dev/res_cibersort.rds")

# Example of top regressions
res_cibersort %>%
  unnest(inference) %>%
  arrange(p.value) %>% 
  slice(1:10) %>%
  unnest(cell_type_proportions) %>%
  ggplot(aes(boot::logit(.proportion), log(PFI.time.2), color=estimate)) +
  geom_point() +
  facet_wrap(~interaction(.cell_type, type))

# Correlation between significance and rarity
library(scales)
logit <- trans_new("logit perc", transform = function(x)qlogis(x), inverse = function(x)plogis(x))

res_cibersort %>%
  unnest(inference) %>%
  mutate(mean_prop = map_dbl(cell_type_proportions, ~ .x %>% pull(.proportion) %>% mean)) %>%
  mutate(n_0s = map_dbl(cell_type_proportions, ~ .x %>% mutate(is_0 = (.proportion==0) %>% factor(levels=c(FALSE, TRUE))) %>% count(is_0, .drop=F ) %>% mutate(tot = sum(n)) %>% mutate(frac = n/tot) %>% filter(is_0==TRUE) %>% pull(frac))) %>%
  ggplot(aes((mean_prop + 1e-5), (p.value))) +
  geom_point() + 
  geom_smooth(method = "lm") + coord_trans(y = logit) + scale_x_log10()
