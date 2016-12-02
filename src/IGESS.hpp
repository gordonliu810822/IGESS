//
//  IGESS.hpp
//  IGESSArma
//
//  Created by DaiMingwei on 16/11/8.
//  Copyright © 2016年 daviddai. All rights reserved.
//

#ifndef IGESS_hpp
#define IGESS_hpp
#include <stdio.h>
#include <math.h>
#include "IGESS_aux.hpp"



struct PairCORAUC{
  double cor;
  double auc;
};


IGESSfit* iGess(float* lpfX, vec y, int P,mat* Z = NULL, mat* lpsummaryinfo = NULL, Options* opt = NULL);
PairCORAUC iGessCV(float* lpfX, vec y, int P, mat* Z = NULL,  mat* lpsummaryinfo = NULL, Options* opt = NULL);
arma::Col<uword> cross_valind(arma::uword N, arma::uword nfold);
#endif /* IGESS_hpp */
