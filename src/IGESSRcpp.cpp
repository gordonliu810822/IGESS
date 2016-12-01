#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
#include <R.h>
#include <Rinternals.h>
#include "Rcpp_aux.hpp"
using namespace Rcpp;
using namespace arma;
#include "IGESS_aux.hpp"
#include "IGESS.hpp"
#include "readPlink.hpp"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

//' @title
//' IGESS
//' @description
//' IGESS Model Training
//'
//' @param X  a N by P matrix of the genotype data, N-sample size, P - SNPs number
//'
//' @param y  a N by 1 matrix of the phenotypes
//'
//' @param Z  a N by M(number of covariates) matrix, the covariates of the N individuals, set to NULL by default
//' @param SS  a K by P matrix for summary statistics, set to NULL by default
//' @param opts  a list with fields(max_iter, display_gap), NULL corresponds to list(max_iter=600, display_gap = 60)
//'
//' @return List of model parameters, refer to \code{wrap_fit} function in Rcpp_auc.hpp for the fields
//'
//' @examples
//' ##Working with no summary statistics, no covariates and options
//' data(IGESSDB)
//' fit <- iGess(Model$X,Model$y)
//'
//'
//' #Working with summary statistics(p-values) and no covariates
//' data(IGESSDB)
//' fit <- iGess(Model$X, Model$y, NULL, Model$AIPVal)
//'
//' #Options Specified
//' data(IGESSDB)
//' opts = list(max_iter=600, display_gap = 100)
//' fit <- iGess(X,y,NULL,SS,opts)
//'
//' @details
//' \code{iGess} fits the iGess model. It requires to provide genotype matrix \code{X}
//' and phenotype matrix \code{y}, while users can also provide covariates data to \code{Z},
//' the summary statistics to \code{SS}, and the settings to \code{opts}.  these last three
//' parameters \code{Z}, \code{SS},\code{opts} could be set to NULL.
//' @export
// [[Rcpp::export]]
RcppExport SEXP iGess(arma::fmat& X, arma::vec& y, SEXP Z = R_NilValue, SEXP  SS = R_NilValue, SEXP opts = R_NilValue){
      Options* lp_opt = NULL;
      if(!Rf_isNull(opts)){
         Rcpp::List opt(opts);
         lp_opt = new Options(opt["max_iter"],opt["display_gap"]);
      }
      mat* Z_ = NULL;
      convert2mat(Z_, Z);
      Mat<double>* Summary = NULL;
      convert2mat(Summary, SS);
      if( Summary != NULL ){
        if(Summary -> n_rows > Summary -> n_cols){
          (*Summary) = Summary -> t();
        }
      }
   //   std::vector<float> datX = Rcpp::as<std::vector<float> >(X);
      float* lpfX = X.memptr();
      cout << "Begin IGESS......" << endl;
      int P = X.n_cols;
      IGESSfit* fit = iGess(lpfX, y, P, Z_,Summary,lp_opt);
      return wrap_fit(fit);
}


//' @title
//' IGESSCV
//' @description
//' evaluating the performance of prediction by cross validation
//'
//' @param X  a N by P matrix of the genotype data, N-sample size, P - SNPs number
//'
//' @param y  a N by 1 matrix of the phenotypes
//'
//' @param Z  a N by M(number of covariates) matrix, the covariates of the N individuals, set to NULL by default
//' @param SS  a K by P matrix for summary statistics, set to NULL by default
//' @param opts  a list with fields(max_iter, display_gap), NULL corresponds to list(max_iter=600, display_gap = 60, n_fold = 5)
//'
//' @return list(auc, cor), auc is for the case-control studies, cor is for quantitative trait studies
//'
//' @examples
//' ##Working with no summary statistics, no covariates and options
//' data(IGESSDB)
//' pair <- iGessCV(Model$X,Model$y)
//'
//'
//' #Working with summary statistics(p-values) and no covariates
//' data(IGESSDB)
//' pair <- iGessCV(Model$X, Model$y, NULL, Model$AIPVal)
//'
//' #Options Specified
//' data(IGESSDB)
//' opts = list(max_iter=600, display_gap = 100)
//' pair <- iGessCV(X,y,NULL,SS,opts)
//'
//' @details
//' \code{iGessCV} evaluating the performance of prediction by cross validation.
//' It requires to provide genotype matrix \code{X}
//' and phenotype matrix \code{y}, while users can also provide covariates data to \code{Z},
//' the summary statistics to \code{SS}, and the settings to \code{opts}.  these last three
//' parameters \code{Z}, \code{SS},\code{opts} could be set to NULL.
//' @export
// [[Rcpp::export]]
RcppExport SEXP iGessCV(arma::fmat& X, arma::vec& y,SEXP Z = R_NilValue, SEXP  SS = R_NilValue, SEXP opts = R_NilValue){
  Options* lp_opt = NULL;
  if(!Rf_isNull(opts)){
    Rcpp::List opt(opts);
    lp_opt = new Options(opt["max_iter"],opt["display_gap"],opt["n_fold"]);
  }
  mat* Z_ = NULL;
  convert2mat(Z_, Z);
  Mat<double>* Summary = NULL;
  convert2mat(Summary, SS);
  if( Summary != NULL ){
    if(Summary -> n_rows > Summary -> n_cols){
      (*Summary) = Summary -> t();
    }
  }
  float* lpfX = X.memptr();
  int P = X.n_cols;
  cout << "Begin CV....." << endl;
  PairCORAUC predict = iGessCV(lpfX, y,P, Z_, Summary, lp_opt);
  Rcpp::List ret;
  ret["auc"] = predict.auc;
  ret["cor"] = predict.cor;
  return ret;
}

//' @title
//' iGessPlink
//' @description
//' IGESS Model Training with plink binary file and summary statistics file
//'
//' @param genoplinkfile  the plink file for the genotype data
//'
//' @param summaryfile  the summary statistics file for the p-values
//'
//' @param configfile  the config file for summary statistics to specify the columns of pvalues
//'
//' @return List of model parameters, refer to \code{wrap_fit} function in Rcpp_auc.hpp for the fields
//'
//' @examples
//' # file name and its format
//' # the required binary files are CD.bed, CD.bim, CD.fam, \code{genoplinkfile} should be \code{CD} in its own directory.
//' # for \code{summaryfile}, take the headers below for example
//' CHR SNP POS major_al minor_al wtccc_info narac1_info narac2_info eira_info canada_info brass_info
//' wtccc_z narac1_z narac2_z eira_z canada_z brass_z cases_MM cases_Mm cases_mm controls_MM controls_Mm
//' controls_mm meta_OR OR_95%CI_lo OR_95%CI_up meta_z meta_2tP CochranQ Q_Pval
//'
//' #we should specify the column names for the summary statistics required and the column names for chromsome and the SNP name
//'
//' #the content of \code{configfile} is like following
//' zvalue=narac1_z,canada_z,brass_z
//' snp=SNP
//' chr=CHR
//'
//' #zvalue should be replaced by pvalue if p-values are provided
//'
//'
//'
//' #Working with no summary statistics
//' iGessPlink(genoplinkfile)
//'
//'
//' #Working with summary statistics
//' iGessPlink(genoplinkfile, summaryfile, configfile)
//'
//'
//' @details
//' \code{iGessPlink} fits the iGess model. It requires to provide binary plink file name to \code{genoplinkfile}, the binary file
//' should be SNP-major type in the current version of IGESS, while the users can also provide
//' summary statistics file name to \code{summaryfile},
//' the config file \code{configfile} is a must when \code{summaryfile} is provided.
//'
//' \code{iGessPlink} take the overlap of genotype data and summary statistics with respect to the SNPs
//'
//' @export
// [[Rcpp::export]]
RcppExport SEXP iGessPlink(Rcpp::String genoplinkfile, Rcpp::String summaryfile, Rcpp::String configfile){
  GenoInfo obj(genoplinkfile);
  Summary summary(summaryfile, configfile);
  summary.cal_overlap(obj);
  fmat Xf = conv_to<fmat>::from(obj.X);
  if( summary.lpsummary != NULL ){
    summary.lpsummary -> replace(0,1e-12);
    if(summary.lpsummary -> n_rows > summary.lpsummary -> n_cols){
      (*summary.lpsummary) = summary.lpsummary -> t();
    }
    cout <<"Number of SNPs = " << summary.lpsummary -> n_cols << endl;
    cout <<"Number of GWAS = " << summary.lpsummary -> n_rows << endl;
  }
  obj.X.clear();
  int P = Xf.n_cols;
  float* lpXf = Xf.memptr();
  IGESSfit* fit = iGess(lpXf, obj.y, P, NULL, summary.lpsummary,NULL);
  return wrap_fit(fit);
}


//' @title
//' iGessPlinkCV
//' @description
//' evaluating the performance of prediction by cross validation with plink binary file and summary statistics file
//'
//' @param genoplinkfile  the plink file for the genotype data
//'
//' @param summaryfile  the summary statistics file for the p-values
//'
//' @param configfile  the config file for summary statistics to specify the columns of pvalues
//'
//' @return list(auc, cor), auc is for the case-control studies, cor is for quantitative trait studies
//' @examples
//' # file name and its format
//' # the required binary files are CD.bed, CD.bim, CD.fam, \code{genoplinkfile} should be \code{CD} in its own directory.
//' # for \code{summaryfile}, take the headers below for example
//' CHR SNP POS major_al minor_al wtccc_info narac1_info narac2_info eira_info canada_info brass_info
//' wtccc_z narac1_z narac2_z eira_z canada_z brass_z cases_MM cases_Mm cases_mm controls_MM controls_Mm
//' controls_mm meta_OR OR_95%CI_lo OR_95%CI_up meta_z meta_2tP CochranQ Q_Pval
//'
//' #we should specify the column names for the summary statistics required and the column names for chromsome and the SNP name
//'
//' #the content of \code{configfile} is like following
//' zvalue=narac1_z,canada_z,brass_z
//' snp=SNP
//' chr=CHR
//'
//' #zvalue should be replaced by pvalue if p-values are provided
//'
//'
//'
//' #Working with no summary statistics
//' iGessPlink(genoplinkfile)
//'
//'
//' #Working with summary statistics
//' iGessPlink(genoplinkfile, summaryfile, configfile)
//'
//'
//' @details
//' \code{iGessPlinkCV} evaluating the performance of prediction by cross validation. It requires to provide binary plink file name to \code{genoplinkfile}, the binary file
//' should be SNP-major type in the current version of IGESS, while the users can also provide
//' summary statistics file name to \code{summaryfile},
//' the config file \code{configfile} is a must when \code{summaryfile} is provided.
//'
//' \code{iGessPlinkCV} take the overlap of genotype data and summary statistics with respect to the SNPs
//'
//' @export
// [[Rcpp::export]]
RcppExport SEXP iGessPlinkCV(Rcpp::String genoplinkfile, Rcpp::String summaryfile, Rcpp::String configfile){
  GenoInfo obj(genoplinkfile);
  Summary summary(summaryfile, configfile);
  summary.cal_overlap(obj);
  fmat Xf = conv_to<fmat>::from(obj.X);
  cout <<"N"<<Xf.n_rows <<" P="<<Xf.n_cols << endl;
  obj.X.clear();
  if( summary.lpsummary != NULL ){
    summary.lpsummary -> replace(0,1e-12);
    cout <<"max="<<max(max(*summary.lpsummary)) << endl;
    cout <<"min="<<min(min(*summary.lpsummary)) << endl;
    if(summary.lpsummary -> n_rows > summary.lpsummary -> n_cols){
      (*summary.lpsummary) = summary.lpsummary -> t();
    }
    cout <<"Number of SNPs = " << summary.lpsummary -> n_cols << endl;
    cout <<"Number of GWAS = " << summary.lpsummary -> n_rows << endl;
  }
  int P = Xf.n_cols;
  float* lpXf = Xf.memptr();

  PairCORAUC predict = iGessCV(lpXf, obj.y, P, NULL, summary.lpsummary,NULL);
  Rcpp::List ret;
  ret["auc"] = predict.auc;
  ret["cor"] = predict.cor;
  return ret;
}

//' @title
//' read_data
//' @description
//' read the binary plink file or correpsonding summary statistics file into objects
//'
//' @param genoplinkfile  the plink file for the genotype data
//'
//' @param summaryfile  the summary statistics file for the p-values
//'
//' @param configfile  the config file for summary statistics to specify the columns of pvalues
//'
//' @return list(X, y) - genotype data and its phenotypes or list(X,y,pvalues) - X, y and its corresponding pvalues
//' @examples
//' # file name and its format
//' # the required binary files are CD.bed, CD.bim, CD.fam, \code{genoplinkfile} should be \code{CD} in its own directory.
//' # for \code{summaryfile}, take the headers below for example
//' CHR SNP POS major_al minor_al wtccc_info narac1_info narac2_info eira_info canada_info brass_info
//' wtccc_z narac1_z narac2_z eira_z canada_z brass_z cases_MM cases_Mm cases_mm controls_MM controls_Mm
//' controls_mm meta_OR OR_95%CI_lo OR_95%CI_up meta_z meta_2tP CochranQ Q_Pval
//'
//' #we should specify the column names for the summary statistics required and the column names for chromsome and the SNP name
//'
//' #the content of \code{configfile} is like following
//' zvalue=narac1_z,canada_z,brass_z
//' snp=SNP
//' chr=CHR
//'
//' #zvalue should be replaced by pvalue if p-values are provided
//'
//' #Working with no summary statistics
//' read_data(genoplinkfile)
//'
//'
//' #Working with summary statistics
//' read_data(genoplinkfile, summaryfile, configfile)
//'
//'
//' @details
//' \code{read_data} evaluating the performance of prediction by cross validation. It requires to provide binary plink file name to \code{genoplinkfile}, the binary file
//' should be SNP-major type in the current version of IGESS, while the users can also provide
//' summary statistics file name to \code{summaryfile},
//' the config file \code{configfile} is a must when \code{summaryfile} is provided.
//'
//' \code{read_data} take the overlap of genotype data and summary statistics with respect to the SNPs
//'
//' @export
// [[Rcpp::export]]
RcppExport SEXP read_data(Rcpp::String genoplinkfile,  SEXP  summaryfile_ = R_NilValue , SEXP  configfile_ = R_NilValue){
  Rcpp::String* summaryfile = NULL;
  List ret;
  if(!Rf_isNull(summaryfile_)){
    summaryfile = new Rcpp::String(summaryfile_);
    Rcpp::String configfile(configfile_);
    std::ifstream ifs(((string)configfile).c_str());
    if(!ifs.is_open()){
      cout <<"Config file is not properly provided!" << endl;
      ifs.close();
      return ret;
    }
    ifs.close();
  }
  GenoInfo obj(genoplinkfile);
  if(summaryfile != NULL){
    Rcpp::String configfile(configfile_);
    Summary summary(*summaryfile, configfile);
    summary.cal_overlap(obj);
    if( summary.lpsummary != NULL ){
      summary.lpsummary -> replace(0,1e-12);
      if(summary.lpsummary -> n_rows > summary.lpsummary -> n_cols){
        (*summary.lpsummary) = summary.lpsummary -> t();
      }
      cout <<"Number of Snps = " << summary.lpsummary -> n_cols << endl;
      cout <<"Number of GWAS = " << summary.lpsummary -> n_rows << endl;
      ret["pvalues"] = *summary.lpsummary;
    }
  }
  fmat Xf = conv_to<fmat>::from(obj.X);
  ret["X"] = Xf;
  ret["y"] = obj.y;
  return ret;
}

//' @title
//' iGessPredict
//' @description
//' predict the output of the given X and Covariates
//'
//' @param fit  List of the IGESSfit fields
//'
//' @param X  the given data matrix
//'
//' @param Z  covariates, set to be NULL by default
//'
//' @return y - predict of the phenotypes
//' @examples
//'
//'
//' @details
//' \code{iGessPredict} predict the output of the given X and Covariates
//'
//' @export
// [[Rcpp::export]]
RcppExport SEXP iGessPredict(SEXP fit_,  arma::mat& X, SEXP  Z = R_NilValue){
  vec ypred(X.n_rows);
  IGESSfit* fit = NULL;
  if(!Rf_isNull(fit_)){
    Rcpp::List fitList(fit_);
    fit = new IGESSfit(fitList["gammas"], fitList["mu"], fitList["cov"]);
    mat* Z_ = NULL;
    convert2mat(Z_, Z);
    ypred = fit -> predict(&X, Z_);
  }else{
    cout << "Invalid input of iGESS fit!" << endl;
    return Rcpp::wrap(ypred);
  }
  return Rcpp::wrap(ypred);
}

// [[Rcpp::export]]
double calaucRcpp(arma::vec label, arma::vec pred){
    return calaucRcpp(label, pred);
}

