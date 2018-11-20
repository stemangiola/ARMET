![ARMET](https://github.com/stemangiola/ARMET/blob/master/armet_logo.png?raw=true)

ARMET-tc infers rates of changes in tissue composition acros a covariate of interest (e.g., treatment status, time or disease grade). 

Input (see below):

- Data frame of tissue gene counts
- Design matrix
- Covariate of interest

```R

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

Output:

- Matrix of proportions of cell type for each provided sample
- Coefficients of change (slope)

# Mode of usage:

ARMET-tc is tipically used with a design matrix as you would use a linear model of differential expression estimator (e.g., Limma or edgeR)

ARMET-tc can be also uses without any design matrix, in this case the only output would be the matrix of proportions of cell type for each provided sample. 

IMPORTANT: in the latter case the algorithm assumes that the query experiment includes homogeneous samples. That is, coming from similar tissue. For example, querying samples of blood together is a good choice, whereas querying samples of pure macrophages, together with samples of pure t-cells is not a good design. In the latter case those two groups should be queried separately.

# Usage:

```R
# Installation
# For a fresh Ubuntu system

###############################################
# Install R/3.5.1
###############################################

sudo apt install apt-transport-https software-properties-common
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
sudo apt update
sudo apt-get install r-base
sudo apt-get install libssl-dev libcurl4-openssl-dev libxml2-dev

# For an initialised machine
###############################################
# Dependencies
###############################################
For non C++14 native machines

1) create ~/.R/Makevars
2) write in it

CXX14 = g++ # or clang++ if you have that
CXX14FLAGS = -O3 -Wno-ignored-attributes

###############################################
# Dependencies
###############################################
library(devtools)  
source("https://bioconductor.org/biocLite.R") 
biocLite("limma")
biocLite("edgeR")
install.packages("data.tree")
install.packages("abind")
install_github("mjskay/tidybayes")  
install.packages("future")
install.packages("reshape")
install.packages("StanHeaders") 
###############################################
install_github("stemangiola/ARMET", args = "--preclean", build_vignettes = FALSE)  
if("package:ARMET" %in% search()) detach("package:ARMET", unload=TRUE, force=TRUE)
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
ARMET_plotPolar(results)
```

![ARMET](https://github.com/stemangiola/ARMET/blob/master/armet_polar.png?raw=true)

```
# Input data shape
test_data

```

