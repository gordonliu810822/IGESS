//
//  IGESS_aux.cpp
//  IGESSArma
//
//  Created by DaiMingwei on 16/11/8.
//  Copyright © 2016年 daviddai. All rights reserved.
//

#include "IGESS_aux.hpp"

double calauc(arma::vec label, arma::vec pred){
    double auc = 0;
    double m = mean(label);
    vec label2 = label;
    label2(find(label >= m)).fill(1);
    label2(find(label <= m)).fill(0);
    label = label2;
    uword N = pred.size();
    uword N_pos = sum(label);
    uvec  idx = sort_index(pred,"descend");
    vec above = linspace<vec>(1, N, N) - cumsum(label(idx));
    auc = (1 - sum(above % label(idx)) / (N_pos * (N-N_pos)));
    auc = auc > 0.5?auc:(1-auc);
    return auc;
}

void igess_update(float* x_j, double* gamma, double* mu, double d, double s, double* xty_pt, double logPi, double sigma2beta, double sigma2e, int N, double xy, double* ytilde_pt, double* lpsummay, vec* lpparam){

    double r = (*gamma) * (*mu);

    (*mu) = s / sigma2e * (xy + d * r -  dotX(x_j, ytilde_pt, N));

    double SSR_logratio(0);

    SSR_logratio = (*mu) * (*mu) / (2 * s) + 0.5 * log( s / sigma2beta);

    double ai = 0; //additional information provided by the summary statistisc
    if (lpsummay != NULL && lpparam != NULL) {
        uword K = lpparam -> size();
        double* lppara = lpparam -> memptr();
        for (uword i = 0; i < K; i++) {
            ai += log(lppara[i]) + (lppara[i] - 1) * log(lpsummay[i]);
        }
        SSR_logratio += ai;
    }


    (*gamma) = 1/(1+exp(-(logPi + SSR_logratio)));

    double rnew = (*gamma) * (*mu);

    addX (ytilde_pt, rnew - r, x_j, N);

}

/*update the parameters of the gamma distributions*/
void update_betaparam(uword P, uword K, double gamma_sum, double * lpsummary, double* lpgamma, vec* lpparams)
{
    if(K == 0) return;
    vec alphalogpvec(K);
    alphalogpvec.fill(0);
    double* lpalphalogp = alphalogpvec.memptr();
    for(uword k=0; k < K; k++)
    {
        for (uword j = 0; j < P; j++)
        {
            double pvalue = *(lpsummary + K * j + k);
            lpalphalogp[k] += (*(lpgamma + j))*(-log(pvalue));
        }
        (*lpparams)[k] = gamma_sum / lpalphalogp[k];
    }

}


vec IGESSfit::predict(mat* X, mat* Z)
{
    uword N = X -> n_rows;
    mat Z0(N, 1, fill::ones);
    if ( Z == NULL ){
        Z = &Z0;
    }else{
        *Z = join_horiz(Z0, *Z);
    }
    vec yhat = (*Z) * this->cov
    + (*X) * conv_to<vec>::from(this->gammas % this->mu);
    return yhat;
}

void IGESSfit::cal_powerfdr(DataModel* model, double threshold, PerformanceObject* po)
{
    vec gFDR = fdr2FDR(this->gammas);
    uvec ufound = find(gFDR < threshold);
    uword nerr = sum((*model -> labels)(ufound) == 0);
    uword nfound = ufound.size();

    double FDR = (nfound != 0) ? nerr * 1.0 / nfound : 0;

    uword nCausal = sum((*model ->labels) != 0);
    double power =  (ufound.size() - nerr) * 1.0  / nCausal ;

    po -> FDR = FDR;
    po -> power = power;

}

double IGESSfit::cal_auc(DataModel* model){
    vec label = conv_to<vec>::from((*model -> labels));
    double auc = calauc(label, this -> gammas);
    return auc;
}




double lb_pvalue(uword P, uword K, double * lpsummary, double* lpgamma, vec* lpparams){
    double lb = 0;
    if(K == 0) return lb;
    double* lpparam = lpparams -> memptr();
    double pvalue = 0;
    for (uword k = 0; k < K; k++) {
        double param_k = *(lpparam + k);
        for (uword j = 0; j < P; j++){
            pvalue =  *(lpsummary + j * K + k);
            lb += (*(lpgamma + j)) * (log(param_k) + (param_k - 1)*log(pvalue));
        }
    }
    return lb;

}

void update_param(uword N, uword P, Vardist vardist, double sumyytilde, vec diagXTX, double& sigma2e, double& sigma2beta, double& pi_p){
    vec term1 = vardist.gamma % (vardist.sigma2beta + square(vardist.mu));
    double term2 = sum((term1 - square(vardist.gamma % vardist.mu)) % diagXTX);

    sigma2e = (sumyytilde + term2) / N;
    double sum_vardist_gamma = sum(vardist.gamma);
    pi_p = sum_vardist_gamma / P;
    sigma2beta = sum(term1) / sum_vardist_gamma;
}

double dotXX (double* x, float* y, uword n) {
    double z = 0;
    for (uword i = 0; i < n; i++)
        z += x[i] * y[i];
    return z;
}

mat MatXfloat(mat& Zt, float* lpfX, int P){
    mat ZtX(Zt.n_rows, P);
    uword N = Zt.n_cols;
    uword n_row = Zt.n_rows;
    Zt = Zt.t();
    double* Z = Zt.memptr();
    for(uword i = 0; i < n_row; i++){
        double* z_row_i = Z + i * N;
        for(uword j = 0; j < P; j++){
            float* x_col_j = lpfX + j*N;
            ZtX.at(i, j) = dotXX(z_row_i, x_col_j, N);

        }
    }
    return ZtX;
}

vec vecXfloat(vec& Zt, float* lpfX, int P){
    vec ZtX(P);
    uword N = Zt.size();
    double* Z = Zt.memptr();
    for(int j = 0; j < P; j++){
            float* x_col_j = lpfX + j*N;
            ZtX.at(j) = dotXX(Z, x_col_j, N);

    }
    return ZtX;
}


void MatSub(float* lpfX, mat& m){
    uword N = m.n_rows;
    uword P = m.n_cols;
    double* m_ptr = m.memptr();
    for(uword i = 0; i < N*P; i++ ){
        *(lpfX + i) -= *(m_ptr + i);
    }
}

//remove the effects of covariates Z
void remove_cov(float* lpfX, int P, vec& y, mat* Z,mat* SZX, mat* SZy)
{
    uword N = y.size();
    mat Z0(N, 1, fill::ones);
    if ( Z == NULL ){
        Z = &Z0;
    }else{
        *Z = join_horiz(Z0, *Z);
    }
    mat Zt = Z -> t();
    mat invZZ = (Zt*(*Z)).i();
    *SZy = invZZ * ( (y.t() * (*Z)).t());
    *SZX = invZZ * MatXfloat(Zt, lpfX, P);//((Zt * (* lpfX)));
    y -= (*Z) * (*SZy);
    mat m = (*Z) * (*SZX);
    MatSub(lpfX, m);

}


//void getDiagsq(double* diag, float* X, Size n, Size p){
//    for (int i = 0; i < p; i++) {
//        float* x = getXColumn(X, i, n);
//        diag[i] = dotXX(x,x,n);
//    }
//}

double dotXX (float* x, float* y, uword n) {
    double z = 0;
    for (uword i = 0; i < n; i++)
        z += x[i] * y[i];
    return z;
}




vec getDiag(float* X, uword P, uword N){
   // uword p = X -> n_cols;
    vec diag(P);
    for (uword i=0; i < P; i++) {
        float* v_i = X + N * i;
        diag[i] = dotXX(v_i, v_i,N);
    }
    return diag;
}


//lower bound for the linear part
double lb_linear(vec ytilde, vec diagXTX, vec y, double sigma2e, Vardist vardist){
    uword n = y.n_elem;
    double lb1 = - (0.5*n)*log(2 * M_PI * sigma2e);
    double lb2 = - 0.5 * sum(square(y-ytilde))/sigma2e;
    double lb3 =
    -0.5*sum( (vardist.gamma % (vardist.sigma2beta + square(vardist.mu)) -
               square(vardist.gamma % vardist.mu)) % diagXTX )/sigma2e;
    return lb1 + lb2 + lb3;

};

vec logpexp(vec x) {
    vec y = x;
    uvec idx = (x < 16);
    y(idx) = log(1 + exp(x(idx)));
    return y;
}

double logpexp(double x) {
    double y = log(1 + exp(x));
    return y;
}

//lower bound for gamma part
double lb_gamma(vec gamma, double log_pi){
    return sum((gamma - 1) * log_pi + (-logpexp(-log_pi)));
};

//lower bound for the klbeta part
double lb_klbeta(Vardist vardist, double sigma2beta){
    double lb = -sum(vardist.gamma % log(vardist.gamma+(vardist.gamma==0)) + (1-vardist.gamma) % log(1-vardist.gamma+(vardist.gamma == 1))) + 0.5*sum(vardist.gamma % (1+log(vardist.sigma2beta / sigma2beta)- (square(vardist.mu) + vardist.sigma2beta)/ sigma2beta ));
    return lb;
}


template <typename T>
double dotX (T* x, double* y, int n) {
    double sum = 0;
    //  #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < n; i++)
        sum += x[i] * y[i];
    return sum;
}


void addX (double* y, double a, float* x, int n) {
    //   #pragma omp parallel for num_threads(4)
    for (int i = 0; i < n; i++)
        y[i] += a * x[i];
}

//convert local fdr to Global FDR
vec fdr2FDR(vec fdr){
    uword M = fdr.size();
    uvec indices = sort_index(fdr);
    vec sort_fdr = fdr(indices);
    vec FDR = cumsum(sort_fdr) / linspace<vec>(1, M, M);
    FDR.elem(find(FDR  > 1)).ones();
    FDR(indices) = FDR;
    return FDR;
}




