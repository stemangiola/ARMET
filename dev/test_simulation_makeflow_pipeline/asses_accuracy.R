library(tidyverse)
library(magrittr)
library(ARMET)
library(tidybulk)

args = commandArgs(trailingOnly=TRUE)


# foreignProp="0.5"
# run="1"
# S = "60"
# slope = 0.05
# whichChanging = "1"
# method ="ARMET"
# input_file ="dev/test_simulation//input__slope_0.05__foreignProp_0.5__S_60__whichChanging_1__run_1.rds"
# file = "dev/test_simulation/output__slope_0.05__foreignProp_0.5__S_60__whichChanging_1__run_1__method_ARMET.rds"
# output_file = "dev/test_simulation/asses__slope_0.05__foreignProp_0.5__S_60__whichChanging_1__run_1__method_ARMET.rds"



foreignProp = args[1]
run = args[2]
S = args[3]
slope = as.numeric(args[4])
whichChanging = args[5]
method = args[6]
input_file = args[7]
file = args[8]
output_file = args[9]

# Pre-load ARMET in case
if(method == "ARMET") ARMET_res = readRDS(file)

decide_if_fp_tp = function(.data, CI, slope){
  .data %>%
    
    # If NA set to 0 and 1
    mutate(
      estimate = if_else(estimate %>% is.na, 0, estimate),
      p.value = if_else(p.value %>% is.na, 1, p.value)
    ) %>%
    
    # If non significant of significant with opposite direction
    mutate(
      fp = alpha_2 == 0 &
        p.value <= CI |
        (alpha_2 != 0 &	p.value < CI & (estimate*alpha_2)<0)
    ) %>%
    
    # If significant but the same direction 
    mutate(tp = alpha_2 != 0 &	p.value < CI & (estimate*alpha_2)>0)  %>%
    
    # # Resolve NA presence that cause error
    # mutate(
    #   fp = if_else(is.na(fp), FALSE, fp),
    #   tp = if_else(is.na(tp), FALSE, tp),
    #   alpha_2 = if_else(is.na(alpha_2), 0, alpha_2)
    # ) %>%
    
    # Filter out accidental fp because of simplex, when cell is really 0 slope
    # Keep the NA because I need for ARMET lower level retrieval
    filter(
      is.na(alpha_2) | !(
        fp & 
        (slope * estimate )<0 &
        alpha_2 == 0
      ) 
    )
}

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
    
    # Apparently the estimate have opposite sign!!!
    mutate(estimate = -estimate) %>%
    
    decide_if_fp_tp(CI, slope)  %>%
    
    # Filter cells that are not present in reference
    filter(alpha_2 %>% is.na %>% `!`)

    
}

check_lower_levels_for_ARMET = function(.data, input_df, CI, slope){
  
  ancestors = 
    ARMET::tree %>% 
    ToDataFrameTypeColFull(TRUE) %>% 
    # LEVEL DEPENDENT
    filter(level_5 %in% (input_df %>% filter(alpha_2!=0))) %>% 
    distinct(level_1, level_2, level_3, level_4) %>% 
    as.data.frame  %>%
    as.character
  
  second_best = 
    .data %>%
    filter(`Cell type category` %in% ancestors ) %>% 
    mutate(divisor = 1:n()) %>%
    
    # Decide what is the penalty
    mutate(tp = divisor/4) %>%
    filter( p.value < CI & (estimate*slope)>0) %>%
    arrange(divisor) %>% 
    dplyr::slice(1) %>%
    pull(tp)
  
  .data %>%
    mutate(tp = ifelse(
      !is.na(alpha_2) &
        alpha_2 != 0 & 
        !tp & 
        length(second_best)>0,
      second_best,
      as.integer(tp)
    )) 
  
}

process_ARMET = function(CI, file, input_file, slope){
 
  # Load input data frame that will be used more than once
  input_df = 
    readRDS(input_file) %>% 
    attr("proportions") %>% 
    distinct(`Cell type category`, alpha_2)
  
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
    left_join(input_df, by = c("Cell type category"  )) %>%
    
    mutate(
      estimate = .value_2,
      p.value = (1-abs(prob_non_0_2))
    ) %>%
    
    #adjust_prob_for_non_hierarchycal_assumption() %>%
    
    decide_if_fp_tp(CI, slope) %>%
    
    # If I don't get the answer, check the other levels
    check_lower_levels_for_ARMET(input_df, CI, slope) %>%
    
    # Filter cells that are not present in reference
    filter(alpha_2 %>% is.na %>% `!`)
  
    # # Calculate
    # mutate(fp = (alpha_2 == 0) &	abs(prob_non_0_2) > CI) %>%
    # mutate(tp = (alpha_2 != 0) &	abs(prob_non_0_2) > CI) %>%
    # # Filter out accidental fp because of simplex
    # filter(!(fp & (slope * .value_2 )<0))  
  
}

adjust_prob_for_non_hierarchycal_assumption = function(.data){
  
  
  adjusted_p = 
    ARMET::tree %>% ToDataFrameTypeColFull() %>%
    select(1:4) %>%
    distinct() %>%
    
    # Attach data
    left_join(
      .data %>% 
        select(`Cell type category`, p.value2 = p.value, estimate2 = estimate) %>%
        distinct, 
      by=c("level_2" = "Cell type category")
    ) %>%
    left_join(
      .data %>% 
        select(`Cell type category`, p.value3 = p.value, estimate3 = estimate) %>%
        distinct, 
      by=c("level_3" = "Cell type category")
    ) %>%
    left_join(
      .data %>% 
        select(`Cell type category`, p.value4 = p.value, estimate4 = estimate) %>%
        distinct, 
      by=c("level_4" = "Cell type category")
    ) %>%
    
    # Correct p.values
    mutate(
      p.value3 = if_else(estimate3 * estimate4 < 0, 1/p.value3, p.value3),
      p.value2 = if_else(estimate2 * estimate4 < 0, 1/p.value2, p.value2)
    ) %>%
    
    mutate(
      p.value = case_when(
        level_2 == level_3 ~ p.value2,
        level_3 == level_4 ~ p.value3 * p.value2,
        TRUE ~ p.value4 * p.value3 * p.value2
    )) %>%
    
    rowwise() %>%
    mutate(p.value = min(1, p.value)) %>%
    
    select(level_4, p.value)
  
  .data %>%
    dplyr::rename(old_p.value = p.value) %>%
    left_join(adjusted_p, by=c("Cell type category" = "level_4"))
}

tibble(
  foreignProp = foreignProp,
  run = run,
  S = S,
  slope = slope,
  whichChanging = whichChanging,
  method = method,
  input_file = input_file,
  file = file
) %>%
  
  mutate(   CI =  list(rev(1-c( seq(0.0, 0.85, 0.05), seq(0.85, 0.9999, 0.001),seq(0.9999, 0.99999999, 0.000001) )))  ) %>%
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






