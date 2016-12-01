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
#include "IGESS_aux.hpp"
#include <RcppArmadillo.h>
#include <math.h>
using namespace arma;

struct PairCORAUC{
  double cor;
  double auc;
};


IGESSfit* iGess(fmat* lpfX, fvec y,fmat* Z = NULL, mat* lpsummaryinfo = NULL, Options* opt = NULL);

//IGESSfit* iGess(fmat* lpfX, fvec y,mat* Z,  mat* lpsummaryinfo = NULL, Options* opt = NULL);


PairCORAUC iGessCV(fmat* lpfX, fvec y, fmat* Z = NULL,  mat* lpsummaryinfo = NULL, Options* opt = NULL);


Col<uword> cross_valind(uword N, uword nfold);
#endif /* IGESS_hpp */
