library(tidyverse)
library(magrittr)
library(ARMET)
library(tidybulk)

args = commandArgs(trailingOnly=TRUE)


foreignProp = args[1]
run = args[2]
S = args[3]
slope = args[4]
whichChanging = args[5]
method = args[6]
input_file = args[7]
file = args[8]
output_file = args[9]

# Pre-load ARMET in case
if(method == "ARMET") ARMET_res = readRDS(file)

process_third_party = function(CI, file, input_file, slope){
  readRDS(file) %>% 
    
    # Add truth
    dplyr::select(.cell_type, estimate, 	p.value) %>%
    left_join(
      readRDS(input_file) %>% 
        attr("proportions") %>% 
        distinct(`Cell type category`, alpha_2),
      by = c(".cell_type" = "Cell type category"  )
    ) %>%
    drop_na  %>%
    
    # Calculate
    mutate(fp = alpha_2 == 0 &	p.value < CI) %>%
    mutate(tp = alpha_2 != 0 &	p.value < CI)  %>%
    
    # Filter out accidental fp because of simplex
    filter(!(fp & (slope * estimate )<0)) 
}

process_ARMET = function(CI, file, input_file, slope){
  #readRDS(file) %$% ## NOW I LOAD THIS BEFORE HAND ONLY ONCE AS IT WAS SLOWING DOWN EVERYTHING
  ARMET_res %$%
    proportions %>% 
    #get_CI(CI) %>%
    
    mutate(.value_2 = map_dbl(
      draws_cens,
      ~ .x %>% 
        filter(A == 2) %>%
        pull(.value ) %>%
        mean
    )) %>%
    
    mutate(prob_non_0_2 = map_dbl(
      draws_cens,
      ~ .x %>% 
        filter(A == 2) %>% 
        ARMET:::draws_to_prob_non_zero()
    )) %>%
    
    # Add truth
    dplyr::select(`Cell type category`, .value_2, prob_non_0_2) %>%
    left_join(
      readRDS(input_file) %>% 
        attr("proportions") %>% 
        distinct(`Cell type category`, alpha_2),
      by = c("Cell type category"  )
    ) %>%
    drop_na  %>%
    
    # Calculate
    mutate(fp = (alpha_2 == 0) &	abs(prob_non_0_2) > CI) %>%
    mutate(tp = (alpha_2 != 0) &	abs(prob_non_0_2) > CI) %>%
    
    # Filter out accidental fp because of simplex
    filter(!(fp & (slope * .value_2 )<0))  
  
}

tibble(
  foreignProp = args[1],
  run = args[2],
  S = args[3],
  slope = as.numeric(args[4]),
  whichChanging = args[5],
  method = args[6],
  input_file = args[7],
  file = args[8],
) %>%
  
  # Add CI interval
  mutate(
    CI = 
      if_else(
        method=="ARMET",
        list(c( seq(0.0, 0.85, 0.05), seq(0.85, 0.9999, 0.001))),
        list(rev(1-c( seq(0.0, 0.85, 0.05), seq(0.85, 0.9999, 0.001),seq(0.9999, 0.99999999, 0.000001) )))
      )
  ) %>%
  unnest(CI) %>%
  
  # Calculate
  mutate(fp_tp = pmap(
    list(CI, file, input_file, slope, method),
    ~ ..5 %>% purrr::when(
      (.) == "ARMET" ~ process_ARMET(..1, ..2, ..3, ..4),
      ~ process_third_party(..1, ..2, ..3, ..4)
      
    ))) %>%
  unnest(fp_tp) %>%

  # Save
  saveRDS(output_file, compress = "xz")






