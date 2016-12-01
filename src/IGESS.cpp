//
//  IGESS.cpp
//  IGESSArma
//
//  Created by DaiMingwei on 16/11/8.
//  Copyright © 2016年 daviddai. All rights reserved.
//

#include "IGESS.hpp"
//#include <libiomp/omp.h>


IGESSfit* iGess(float* lpfX, vec y, int P, mat* Z,  mat* lpsummaryinfo, Options* opt ){
    cout << "Start of the main function...." << endl;
    uword N = y.n_rows;
    //uword P = lpfX -> n_cols;

    mat SZX;
    mat SZy;
    remove_cov(lpfX, P, y, Z, &SZX, &SZy);

    uword K = lpsummaryinfo != NULL ? (lpsummaryinfo -> n_rows) : 0;
    cout <<"K=" << K << endl;
    opt = opt != NULL ? opt : new Options();
    int max_iter = opt -> max_iter;
    int display_gap = opt -> display_gap;

    vec xty = vecXfloat(y, lpfX, P);//(y.t() * (*lpfX)).t();
    vec diagXTX = getDiag(lpfX, P, N);
    double pi_p = 0.01; //pi for prior proportion
    double mu0 = 0;
    double alpha0 = 0.5; //initial parameters of beta distribtuion

    Vardist vardist(P, mu0, pi_p);

    double sigma2e = var(y) / 2;
    double sigma2beta = sigma2e;

    vec beta = vardist.gamma % vardist.mu;
    vec ytilde = vec(N);//(*lpfX) * beta;
    ytilde.zeros();

    Col<double>* lpparams = NULL;
    if ( K > 0 ){
        lpparams = new Col<double>(K);
        lpparams -> fill(alpha0); //parameters for beta distribtuions
    }

    double L0 = -INFINITY;
    double L = 0;
    double* lpgamma = vardist.gamma.memptr();
    double* lpmu = vardist.mu.memptr();
    float* lpX = lpfX;
    double* lpd = diagXTX.memptr();
    double* lpytilde = ytilde.memptr();
    double* lpxy = xty.memptr();
    double* lpsummary = lpsummaryinfo != NULL ? lpsummaryinfo -> memptr() : NULL;
    uword iter = 0;

    for (iter = 0; iter < max_iter; iter ++) {
        clock_t t1 = clock();
        if(iter == 0)  cout <<"Begin Iterations" << endl;
        double logPi = log(pi_p / (1 - pi_p));
        double sigma2e_Sigma2beta = sigma2e / sigma2beta;
        vec xxsigma = diagXTX + sigma2e_Sigma2beta;
        vardist.sigma2beta = sigma2e / xxsigma;
        double* S = vardist.sigma2beta.memptr();

        double gamma_sum = 0;
        for (int j = 0; j < P; j++) {
            igess_update(lpX + N*j, lpgamma + j, lpmu + j, *(lpd + j), *(S+j), lpxy, logPi, sigma2beta, sigma2e, (int)N,  *(lpxy + j), lpytilde, lpsummary + K * j, lpparams);
            gamma_sum += *(lpgamma + j);

        }

        update_betaparam(P, K, gamma_sum, lpsummary, lpgamma, lpparams);

        update_param(N, P, vardist, sum(square(y-ytilde)), diagXTX, sigma2e, sigma2beta, pi_p);


        L = lb_linear(ytilde, diagXTX,y, sigma2e, vardist) + lb_gamma(vardist.gamma, logPi) +lb_klbeta(vardist, sigma2beta) + lb_pvalue(P, K, lpsummary, lpgamma, lpparams);

        if(iter % display_gap == 0){
            printf("%llu iteration,L=%f sigma2e = %f sigma2beta = %f time = %f \n",iter, L, sigma2e, sigma2beta, (clock() - t1)*1.0/CLOCKS_PER_SEC);
        }
        if(L < L0){
            printf("Lowerbound decreasing,Error at iteration %llu th iteration, diff=%g",iter,L-L0);
            break;
        }else if(L - L0 < 1e-5){
            printf("Converge at %lluth iteration",iter);
            break;
        }
        L0 = L;

    }
    cout <<"L=" << L << endl;
    mat cov = SZy - SZX * conv_to<vec>::from(vardist.gamma % vardist.mu);
    IGESSfit* fit = new IGESSfit(N, P,  K, iter, L,  sigma2e, sigma2beta, pi_p, vardist.gamma, vardist.mu
                                 , vardist.sigma2beta, lpparams,  ytilde, cov);
    return fit;
}

/******************************************************************
 Function:       IGESSCV
Description:    Calculate the correlation and auc of the cross validation of the input data
Calls:
Input:
mat* lpfX,   //(*lpfX) is a N by P matrix of float,
N denotes the number of samples,
P is the number of SNPs
vec y,       //y is N by 1 vector of float corresponding to phenotype of each individual
mat* Z,      // (*Z) is a N by M matrix of float,  M is the number of covariates for each individual
//default value : NULL
mat* lpsummaryinfo,   //(*lpsummaryinfo) is a P by K matrix, K is the number of GWAS with Summary Statistics
//default value : NULL
Options* opt          // the options for the functions
//default value : NULL
Return:
 PairCORAUC      // A Struct of the value of AUC and Correlation
******************************************************************/
PairCORAUC iGessCV(float* lpfX, vec y, int P, mat* Z,  mat* lpsummaryinfo, Options* opt){
    opt = opt != NULL ? opt : new Options();
    uword nfold = opt -> n_fold;
    uword N = y . n_rows;
   // uword P = lpfX -> n_cols;
    Col<uword> indices = cross_valind(N, nfold);
    mat SZX;
    mat SZy;
    vec ylabel = y;
    remove_cov(lpfX, P, y, Z, &SZX, &SZy);
    vec predY(N);
    Mat<float> Xf = Mat<float>(lpfX, N, P, false);
    for (uword i = 1; i <= nfold; i++) {
        cout <<"fold i=" << i << endl;
        Col<uword> train_idx = find(indices != i);
        Col<uword> test_idx = find(indices == i);
        Mat<float> trainM = Xf.rows(train_idx);
        vec ytrain = y(train_idx);
        Mat<double> testM = conv_to<mat>::from(Xf . rows(test_idx));
        vec ytest = y(test_idx);
        float* trainX = trainM.memptr();
        IGESSfit* f = iGess(trainX, ytrain, P,  NULL, lpsummaryinfo, opt);
        vec predy = f -> predict(&testM);
        predY.elem(test_idx) = predy;
        double accuracy = calauc(conv_to<vec>::from(ylabel(test_idx)), conv_to<vec>::from(predy));
        cout << "accuracy = " << accuracy << endl;
        delete f;

    }
    double accuracy = calauc(conv_to<vec>::from(ylabel), conv_to<vec>::from(predY));
    double corr = as_scalar(cor(y, predY));
    PairCORAUC pair;
    pair.auc = accuracy;
    pair.cor = corr;
    return pair;
}


Col<uword> cross_valind(uword N, uword nfold){
    arma_rng::set_seed_random();
    Col<uword> vec_n = shuffle(linspace < Col <uword> >(1, N, N));
    Col<uword> indices(N);
    indices.fill(nfold);
    for(uword n = 1; n <= nfold-1; n++){
        Col<uword> in = vec_n.rows((n-1)*N/nfold,n*N/nfold - 1);
        indices.elem(in - 1 ).fill(n);
    }
    return indices;
}
