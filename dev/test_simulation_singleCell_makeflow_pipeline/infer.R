library(EPIC)
library(tidyverse)
library(magrittr)
library(ARMET)
library(tidybulk)



args = commandArgs(trailingOnly=TRUE)
input_file = args[1]
output_file = args[2]
method = args[3]

get_survival_X = function(S){
  readRDS("dev/PFI_all_cancers.rds") %>%
    filter(PFI.2 == 1 & !is.na(PFI.time.2) & PFI.time.2 > 0) %>%
    select(real_days = PFI.time.2 ) %>%
    mutate(real_days = real_days %>% scale(center = F) %>% as.numeric) %>%
    sample_n(S) %>%
    mutate(sample = sprintf("S%s", 1:n())) %>%
    mutate(alive = sample(0:1, n(), replace = T)) %>%
    mutate(days = ifelse(alive==1, real_days/2, real_days) ) %>%
    mutate(intercept = 1)
}

get_noiseless_harmonised = function(){
  
  mix_base_unharmonized = readRDS("dev/mix_base_noiseless.RDS")
  
  my_markers =
    ARMET::ARMET_ref %>%
    
    left_join(ARMET::n_markers, by = c("ct1", "ct2")) %>%
    filter_reference(
      mix_base_unharmonized %>%
        filter(level == 3) %>%
        distinct(`Cell type category`, symbol, `count normalised bayes`) ,
      ARMET::n_markers
    ) %>% distinct(level, symbol)
  
  # level 1
  abundance_1 =
    my_markers %>% filter(level == 1) %>%
    left_join(mix_base_unharmonized) %>%
    select(level_2, symbol,  `count normalised bayes 1` =`count normalised bayes`)
  
  abundance_2 =
    my_markers %>% filter(level == 2) %>%
    left_join(mix_base_unharmonized) %>%
    select(level_3, symbol,  `count normalised bayes 2` =`count normalised bayes`)
  
  # Now this is noiseless for the ancestor markers so also for ARMET that rely on hierarchy
  mix_base_unharmonized %>%
    filter(level==3) %>%
    left_join(abundance_2) %>%
    left_join(abundance_1) %>%
    mutate(`count normalised bayes 2` = ifelse(`count normalised bayes 1` %>% is.na, `count normalised bayes 2`, `count normalised bayes 1`)) %>%
    mutate(`count normalised bayes` = ifelse(`count normalised bayes 2` %>% is.na, `count normalised bayes`, `count normalised bayes 2`)) %>%
    select(level_2, level_3, level_4, `Cell type category`, level, sample, symbol, `count normalised bayes`, `house keeping`)
  
}

# get_noiseless_harmonised() %>%
#   distinct(`Cell type category`, symbol, `count normalised bayes`) %>%
#   spread(`Cell type category`, `count normalised bayes`) %>%
#   nanny::as_matrix(rownames = symbol) %>%
#   as.data.frame %>% saveRDS("dev/test_simulation_makeflow_pipeline/third_party_reference.rds", compress = "xz")

readRDS(input_file) %>%
  
  when(
    
    # ARMET
    method == "ARMET" ~ {
      mix = (.)
      
      result =
        mix %>%
        ARMET_tc(
          ~ censored(days, alive),
          sample, symbol, `count mix`,
          prior_survival_time =  
            mix %>% 
            distinct(sample) %>% 
            nrow %>%
            get_survival_X() %$% 
            real_days %>% 
            as.numeric, 
          iterations = 500, 
          sampling_iterations = 300
        )  %>%
        ARMET_tc_continue(2) %>%
        ARMET_tc_continue(3) %>%
        ARMET_tc_continue(4)
      
      # library(furrr)
      # plan(multisession, workers=4)
      
      result$proportions = result %>% add_cox_test(relative = FALSE)
      
      # Make object lighter
      result_light = result$proportions %>% select(-draws, -rng_prop, -rng_mu) 
      list(proportions = result_light, fit = result$internals$fit)
    },
    
    # If any other method
    ~ (.) %>%
      mutate(dead = !alive) %>%
      tidybulk::test_differential_cellularity(
        survival::Surv(days, dead) ~ .,
        sample, symbol, `count mix`,
        reference = readRDS("dev/test_simulation_makeflow_pipeline/third_party_reference.rds"),
        method = method
      ) 
  ) %>%
  
  # Save
  saveRDS(output_file, compress = "xz")






