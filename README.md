![alt tag](https://github.com/stemangiola/ARMET/blob/master/armet_logo.png?raw=true)

Usage:

install.packages("devtools")  
library(devtools)  
options(buildtools.check = function(action) TRUE)  
install_github("stemangiola/ARMET", args = "--preclean", build_vignettes = FALSE, auth_token = "37c5c6238136a6804d336d9a7078eece993ce870", password="x-oauth-basic")  
library(ARMET)  
tissue_composition = ARMET_tc(test_mix)  
boxplot(tissue_composition$proportions)


