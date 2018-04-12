![alt tag](https://github.com/stemangiola/ARMET/blob/master/armet_logo.png?raw=true)

Usage:

# Installation
library(devtools)  
install_github("stemangiola/ARMET", args = "--preclean", build_vignettes = FALSE, auth_token = "37c5c6238136a6804d336d9a7078eece993ce870", password="x-oauth-basic")  
library(ARMET) 

# Inference
results = 
  ARMET_tc(
    mix =          test_data$mix, 
    my_design =    test_data$design_matrix, 
    cov_to_test =  test_data$cov_to_test, 
    do_debug = F
   )
   
# Report
ARMET_getFit(results)
ARMET_plotFit(results, "immune_cell")

# Input data shape
test_data

```
$mix
# A tibble: 24,507 x 16
   gene     `TCGA-CC-A8HT` `TCGA-DD-AACZ` `TCGA-DD-A73G` `TCGA-DD-A3A5`
   <fct>             <dbl>          <dbl>          <dbl>          <dbl>
 1 TSPAN6             2984           8280        3508           2783
 2 TNMD                  0              0           1.00           3.00
 3 DPM1               1336           1751         521           1548
 4 SCYL3               565            650         307            301
 5 C1orf112            550            348          73.0          186
 6 FGR                 538            830          62.0           90.0
 7 CFH               15926          38744       70886          42886
 8 FUCA2              7318           9106        4212           2153
 9 GCLC               3799           5108         582           2776
10 NFYA               1593           1633         307            333
# ... with 24,497 more rows, and 11 more variables: `TCGA-DD-A4ND` <dbl>,
#   `TCGA-DD-AACX` <dbl>, `TCGA-DD-AAVY` <dbl>, `TCGA-DD-AAVV` <dbl>,
#   `TCGA-DD-AAVR` <dbl>, `TCGA-DD-AAW1` <dbl>, `TCGA-DD-AAVZ` <dbl>,
#   `TCGA-DD-A4NV` <dbl>, `TCGA-UB-A7MB` <dbl>, `TCGA-DD-A4NK` <dbl>,
#   `TCGA-DD-AAD6` <dbl>

$design_matrix
# A tibble: 15 x 4
   sample       `(Intercept)` relapse   age
   <fct>                <dbl>   <dbl> <dbl>
 1 TCGA-CC-A8HT          1.00    1.00 27334
 2 TCGA-DD-AACZ          1.00    1.00 23230
 3 TCGA-DD-A73G          1.00    0    26949
 4 TCGA-DD-A3A5          1.00    1.00 24288
 5 TCGA-DD-A4ND          1.00    1.00 20782
 6 TCGA-DD-AACX          1.00    1.00 24329
 7 TCGA-DD-AAVY          1.00    0    20709
 8 TCGA-DD-AAVV          1.00    0    20751
 9 TCGA-DD-AAVR          1.00    0    16157
10 TCGA-DD-AAW1          1.00    0    20288
11 TCGA-DD-AAVZ          1.00    0    14005
12 TCGA-DD-A4NV          1.00    0    22438
13 TCGA-UB-A7MB          1.00    1.00  8951
14 TCGA-DD-A4NK          1.00    1.00 29244
15 TCGA-DD-AAD6          1.00    1.00 24185

$cov_to_test
[1] "relapse"

```

