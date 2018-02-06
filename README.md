![alt tag](https://github.com/stemangiola/ARMET/blob/master/armet_logo.png?raw=true)

Usage:

install.packages("devtools")  
library(devtools)  
options(buildtools.check = function(action) TRUE)  
install_github("stemangiola/ARMET")  
library(ARMET)  
tissue_composition = ARMET_tc(test_mix)  
boxplot(tissue_composition$proportions)


