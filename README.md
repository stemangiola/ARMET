sccomp - Outlier-aware and count-based compositional analysis of
single-cell data
================
Stefano Mangiola

<!-- badges: start -->

[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![R build
status](https://github.com/stemangiola/tidyseurat/workflows/R-CMD-check/badge.svg)](https://github.com/stemangiola/tidyseurat/actions/)
<!-- badges: end -->

# Installation

``` r
devtools::install_github("stemangiola/ARMET")
```

# Usage

``` r
library(ARMET)
```

    ## Loading required package: Rcpp

    ## Warning: replacing previous import 'tidyr::extract' by 'rstan::extract' when
    ## loading 'ARMET'

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
data("test_mixture")
data("no_hierarchy_reference")

 estimates = 
    test_mixture |>
    convoluted_glm(
   ~ factor_of_interest,
   .sample = sample,
   .transcript = symbol,
   .abundance = count,
   reference = no_hierarchy_reference 
  )
```

    ## Warning in setup_convolved_lm_NON_hierarchical(.data, .formula = .formula, :
    ## tidybulk says: the data does not have the same number of transcript per sample.
    ## The data set is not rectangular.

    ## Warning in aggregate_duplicated_transcripts_bulk(.data, .sample = !!.sample, :
    ## tidybulk says: for aggregation, factors and logical columns were converted to
    ## character

    ## Converted to characters

    ## factorfactorlogical

    ## Warning in warning_if_data_is_not_rectangular(.data, !!.sample, !!.transcript, :
    ## tidybulk says: the data does not have the same number of transcript per sample.
    ## The data set is not rectangular.

    ## No group or design set. Assuming all samples belong to one group.

    ## Warning in warning_if_data_is_not_rectangular(.data, !!.sample, !!.transcript, :
    ## tidybulk says: the data does not have the same number of transcript per sample.
    ## The data set is not rectangular.

    ## Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
    ## Running the chains for more iterations may help. See
    ## https://mc-stan.org/misc/warnings.html#bulk-ess

    ## Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
    ## Running the chains for more iterations may help. See
    ## https://mc-stan.org/misc/warnings.html#tail-ess

    ## # A tibble: 0 × 11
    ## # … with 11 variables: par <chr>, mean <dbl>, se_mean <dbl>, sd <dbl>,
    ## #   2.5% <dbl>, 25% <dbl>, 50% <dbl>, 75% <dbl>, 97.5% <dbl>, n_eff <dbl>,
    ## #   Rhat <dbl>

    ## Warning: Expected 5 pieces. Additional pieces discarded in 126 rows [1, 2, 3, 4,
    ## 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].

    ## Warning: Expected 5 pieces. Additional pieces discarded in 126 rows [1, 2, 3, 4,
    ## 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].

    ## Joining, by = c("Q", "sample")

    ## Warning: Expected 5 pieces. Additional pieces discarded in 246 rows [1, 2, 3, 4,
    ## 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].

    ## Warning: Expected 5 pieces. Additional pieces discarded in 42 rows [1, 2, 3, 4,
    ## 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].

    ## Warning in mask$eval_all_mutate(quo): NAs introduced by coercion

    ## Warning: Expected 5 pieces. Additional pieces discarded in 42 rows [1, 2, 3, 4,
    ## 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].

    ## Warning in mask$eval_all_mutate(quo): NAs introduced by coercion

``` r
 estimates
```

    ## # A tibble: 21 × 5
    ##    cell_type `.median_(Inte…` .median_factor_… `.sd_(Intercep…` .sd_factor_of_i…
    ##    <fct>                <dbl>            <dbl>            <dbl>            <dbl>
    ##  1 endothel…           -0.717          -0.248             0.406            0.456
    ##  2 epitheli…           -0.582          -0.313             0.418            0.480
    ##  3 fibrobla…           -0.677          -0.292             0.458            0.486
    ##  4 mast_cell           -0.804          -0.299             0.485            0.444
    ##  5 b_memory             2.04            2.30              0.396            0.397
    ##  6 b_naive              5.93           -0.215             0.230            0.296
    ##  7 eosinoph…            0.284           0.705             0.486            0.504
    ##  8 monocyte            -0.321          -0.0574            0.447            0.482
    ##  9 neutroph…           -0.165           0.130             0.482            0.556
    ## 10 nk_resti…           -0.391          -0.0737            0.431            0.486
    ## # … with 11 more rows
