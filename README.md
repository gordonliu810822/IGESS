#IGESS

IGESS is a statistical approach to integrating individual level genotype data and summary statistics in Genome Wide Association Studies. 'IGESS' R package provides computationally efficient and user friendly interface to fit and evaluate the IGESS model. It accepts both the R-type data  and binary plink files.

##Usage

The following two help pages provide a good start point for the genetic analysis using IGESS package, including the overview of IGESS package and the example command lines: 

library(IGESS)  
help(package="IGESS")  

##Development 
To install the development version of IGESS, it's easiest to use the 'devtools' package. Note that IGESS depends on the 'Rcpp' and 'RcppArmadillo' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.

install.packages("devtools")  
library(devtools)  
install_github("daviddaigithub/IGESS")  

##References

=======
# IGESS
IGESS(Gentic Analysis integrating individual level data and summary statistics)
>>>>>>> c0a44e88fea7a52e78ceb91619fc72079d069f9b
