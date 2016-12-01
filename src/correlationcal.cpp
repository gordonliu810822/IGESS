#include "correlationcal.hpp"
#include <iostream>
#include <armadillo>
using namespace std;
using namespace arma;
#include "time.h"

//to standarize the vector
void vec_sub(vec& v0, double m_value, double squaresum){
 //    v0 = (v0 - m_value)/squaresum;
    v0 -= m_value;
    v0 /= squaresum;
}


arma::fmat calcorr(arma::Mat<unsigned>* X0, arma::vec index, arma::uvec availIndex, int bandwidth){
  uword M = X0 -> n_cols;
  arma::Mat<unsigned>* X = X0;
  arma::Mat<unsigned> subX;
  if(availIndex.size() != 0 && availIndex.size() < M){
    subX = X0->cols(availIndex);
    X = &subX;
  }

  uword size1 = X->n_rows;
  uword size2 = X->n_cols;

  X -> replace(3, 0); //replace the missing value indicated by 3
  vec meanX(size2);
  vec sqrtsum(size2);
  for (int i=0; i < size2; i++) { //calculate the mean of the vector and sqrt sum
    meanX[i] = sum(X -> col(i))*1.0/size1;
    vec v_i = conv_to<vec>::from(X -> col(i));
    v_i -=  meanX[i];
    mat pd = v_i.t() * v_i;
    sqrtsum[i] = sqrt(pd.at(0));
  }

  vec pos(index.size());
  pos[0] = index[0];
  index = cumsum(index);
  for(int i=1; i <= index.size(); i++){
    pos[i] = index[i];
  }

  arma::fmat corr;//to store the resulted correlation matrix
  corr.zeros(size2, 2 * bandwidth + 1);

  int chrom_idx = 0;
  clock_t t1 = clock();
  /*to store the standarized vector in the near time, the max number is set to bidwidth*/
  mat X_tmp(size1,0);
  vec v_i(size1);
  vec v_j(size1);
  bool end_of_chromsome = false;
  for(int i = 0; i < size2; i++){
    uword snp_idx = i;
    if( availIndex.size() > 0 ){
      snp_idx = availIndex[i];
    }
    while(snp_idx >= pos[chrom_idx]){
      chrom_idx++;
      if(chrom_idx > 22){
         end_of_chromsome = true;
         break;
      }
      X_tmp.reset();
    }
    if(end_of_chromsome){
      break;
    }
    if(X_tmp.n_cols == 0){ //x_i is the first element of some chromsome
      v_i = conv_to<vec>::from(X ->col(i));
      if(sqrtsum[i] > 0)
        vec_sub(v_i,  meanX[i], sqrtsum[i]);
      else{
        v_i.zeros();
      }
    }else{
      v_i = X_tmp.col(0); // extract the first element and then remove
      if(X_tmp.n_cols >= 2){
        X_tmp = X_tmp.cols(1,X_tmp.n_cols - 1);
      }else{
        X_tmp.reset();
      }

    }

    if(i % 10000 == 0){
      cout <<i<< "," ;
      cout <<"Elapsed time is " << (clock()-t1)*1.0/CLOCKS_PER_SEC << endl;
    }

    if(i % 100000 == 0) cout << endl;

    corr(i,bandwidth) = 1;
    for(int j = i + 1; j <= i + bandwidth; j++){
      snp_idx = j;
      if( availIndex.size() > 0 ){
        if(j < availIndex.size())
          snp_idx = availIndex[j];
        else
          break;
      }
      if(snp_idx < pos[chrom_idx]) {
        //extract the v_j from the temp container and insert into the container if it is not in the container
        if(j - i > X_tmp.n_cols){
          v_j = conv_to<vec>::from(X ->col(j));
          if(sqrtsum[j] > 0)
            vec_sub(v_j,  meanX[j], sqrtsum[j]);
          else{
            v_j.zeros();
          }
          X_tmp.insert_cols(X_tmp.n_cols, v_j);
        }else{
          v_j = X_tmp.col(j - i - 1);
        }
        mat dotij = v_i.t() * v_j;
        corr(i, bandwidth + j-i) = dotij.at(0);
      }
      else{
        break;
      }

    }
  }

  for(int i=0; i < size2; i++){
    for(int j = 1; j <= bandwidth && j <= i; j++){
      corr(i, bandwidth - j) = corr(i-j, bandwidth + j);
    }
  }
  return corr;


}

/**
default function to calculate the correlation with no selected index
*/
arma::fmat calcorr(arma::Mat<unsigned>* X, arma::vec index,int bandwidth){
  uvec idx;
  idx.reset();
  return calcorr(X, index, idx ,bandwidth);

}
