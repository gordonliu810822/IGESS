IGESS
=======

IGESS is a statistical approach to integrating individual level genotype data and summary statistics in Genome Wide Association Studies. 'IGESS' R package provides computationally efficient and user friendly interface to fit and evaluate the IGESS model. It accepts both the R-type data  and binary plink files.

Usage
=======

The following two help pages provide a good start point for the genetic analysis using IGESS package, including the overview of IGESS package and the example command lines: 

library(IGESS)  
help(package="IGESS")  

Development 
=======
This R package is developed by Mingwei Dai and Can Yang, and maintained by Can Yang <eeyangc@gmail.com>.

Installation
=======
To install the development version of IGESS, it's easiest to use the 'devtools' package. Note that IGESS depends on the 'Rcpp' and 'RcppArmadillo' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.

install.packages("devtools")  
library(devtools)  
install_github("daviddaigithub/IGESS0")  

References
=======
M. Dai, J. Ming, M., Cai, J. Liu, C. Yang, X. Wan, and Z. Xu. IGESS: A statistical approach to integrating individual level genotype data and summary statistics in genome wide association studies. Bioinformatics. 2017 Sep 15;33(18):2882-2889. doi: 10.1093/bioinformatics/btx314

=======
# IGESS
IGESS(Gentic Analysis integrating individual level data and summary statistics)
>>>>>>> c0a44e88fea7a52e78ceb91619fc72079d069f9b
