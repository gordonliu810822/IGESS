#ifndef Rcpp_aux_hpp
#define Rcpp_aux_hpp
#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "IGESS_aux.hpp"
using namespace arma;
using namespace Rcpp;

//convert the data.frame in R to Matrix
NumericMatrix DF2mat(SEXP input) {
  if(!Rf_isMatrix(input)){
    Function asMatrix("as.matrix");
    return asMatrix(input);
  }else{
    Rcpp::NumericMatrix T(input);
    return T;
  }
}


//convert to SEXP to required Mat<double>
void convert2mat(Mat<double>*& obj, SEXP input){
  if(!Rf_isNull(input)){
    Rcpp::NumericMatrix T = DF2mat(input);
    obj = new Mat<double>(T.begin(), T.nrow(), T.ncol(), false);
  }
}

// //convert to SEXP to required Mat<float>
// void convert2fmat(Mat<float>* obj, SEXP input){
//   if(!Rf_isNull(input)){
//     Rcpp::NumericMatrix T = DF2mat(input);
//     arma::mat Td(T.begin(), T.nrow(), T.ncol(), false);
//     (*obj) = conv_to<fmat>::from(Td);
//   }
// }

RcppExport SEXP wrap_fit(IGESSfit* fit){
  Rcpp::List ret;
  ret["sigma2beta"] = fit -> sigma2beta;
  ret["sigma2e"] = fit -> sigma2e;
  ret["gammas"] = fit -> gammas;
  ret["mu"] = fit -> mu;
  ret["S"] = fit -> S;
  ret["pi"] = fit -> Pi;
  ret["P"] = fit ->P;
  ret["fdr"] = 1 - fit -> gammas;
  ret["cov"] = fit -> cov;
  ret["L"] = fit -> L;
  ret["iter"] = fit -> iter;
  if(fit -> pParam != NULL){
    ret["param_beta"] = (*fit -> pParam);
  }
  return ret;
}





#endif /* IGESS_hpp */
