#ifndef correlationcal_hpp
#define correlationcal_hpp
#include <RcppArmadillo.h>
#include <stdio.h>
using namespace arma;
using namespace std;
// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
arma::fmat calcorr(arma::Mat<unsigned>* X, arma::vec index, int badwidth = 200);
arma::fmat calcorr(arma::Mat<unsigned>* X, arma::vec index, arma::uvec availIndex, int badwidth = 200);

//arma::mat calCorrU(arma::umat X, arma::vec index,int badwidth = 200);


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//
#endif /* correlationcal_hpp */
