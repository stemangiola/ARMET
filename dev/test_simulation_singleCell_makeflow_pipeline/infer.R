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
        ARMET_tc_continue(3) 
      
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
      
      # Deconvolution
      tidybulk::deconvolve_cellularity(
        sample, symbol, `count mix`,
        method=method,
        reference = readRDS("dev/test_simulation_singleCell_makeflow_pipeline/third_party_reference_lv3.rds"),
        action="get"
      ) %>%
    
      # Test
      pivot_longer(
        names_prefix = sprintf("%s: ", method),
        cols = starts_with(method),
        names_to = ".cell_type",
        values_to = ".proportion"
      ) %>%
      
      
      # Adapt cell types to single cells
      mutate(
        .cell_type = case_when(
          .cell_type %in% c("b_memory", "b_naive") ~ "b_cell",
          .cell_type %in% c("nk_primed", "nk_resting") ~ "natural_killer",
          TRUE ~ .cell_type
        )
      ) %>%
      nanny::nest_subset(data = -c(sample, .cell_type)) %>% 
      mutate(.proportion = map_dbl(data, ~ sum(.x$.proportion))) %>% 
      select(-data) %>%
      
      
      # Replace 0s
      mutate(min_proportion = min(.proportion[.proportion!=0])) %>%
      mutate(.proportion_0_corrected = if_else(.proportion==0, min_proportion, .proportion)) %>%
      
      # Test survival
      tidyr::nest(cell_type_proportions = -.cell_type) %>%
      mutate(surv_test = map(
        cell_type_proportions,
        ~ {
          if(pull(., .proportion_0_corrected) %>% unique %>% length %>%  `<=` (3)) return(NULL)
          
          # See if regression if censored or not
          .x %>%
            mutate(.proportion_0_corrected = .proportion_0_corrected  %>% boot::logit()) %>%
            survival::coxph(survival::Surv(days, dead) ~ .proportion_0_corrected, .)	%>%
            broom::tidy() %>%
            select(-term)
        }
      )) %>%
      
      unnest(surv_test, keep_empty = TRUE) 
    
  ) %>%
  
  # Save
  saveRDS(output_file, compress = "xz")






