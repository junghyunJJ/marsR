// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;
using namespace arma;

List combnR_marsRR(int n, int k){
  Function combnR_marsR("combnR_marsR");
  List l = combnR_marsR(Named("n", n), Named("k", k));
  return l;
}

NumericVector orderC(NumericMatrix x) {
  mat X(x.begin(), x.nrow(), x.ncol(), false);
  return(as<NumericVector>(wrap(sort_index( abs(X), "descend" ))) );
}

NumericVector computePvalue(NumericMatrix stat){
  NumericVector stat_v = stat( _ , 0 );
  NumericVector abs_max_stat = {-1 * max(abs(stat_v))};
  NumericVector pval = 2 * pnorm(abs_max_stat);
  return(pval);
} 

// [[Rcpp::export]]
NumericMatrix scaleC(NumericMatrix m) {
  NumericMatrix M(m.nrow(), m.ncol());
  for(int i=0; i<M.ncol(); i++){
    NumericVector sel_col = m(_,i);
    double mean_sel_col = mean(sel_col);
    double sd_sel_col = sd(sel_col);
    M(_,i) = (sel_col - mean_sel_col) / sd_sel_col;
  }
  return M;
}


double fracdmvnorm(mat Z, mat mean, mat R, mat diagC, double NCP, int baseValue) {
  mat newR = R + R * diagC  * R;
  mat ZcenterMean = Z - mean;
  // mat res1 = trans(ZcenterMean) * inv(R) * (ZcenterMean);
  // mat res2 = trans(ZcenterMean) * inv(newR) *  (ZcenterMean);
  mat res1 = trans(ZcenterMean) * solve(R, eye(size(R))) * (ZcenterMean);
  mat res2 = trans(ZcenterMean) * solve(newR, eye(size(newR))) *  (ZcenterMean);
  double v1 = res1(0,0)/2-res2(0,0)/2-baseValue/2;
  return(exp(v1)/sqrt(det(newR))* sqrt(det(R)));
}

//[[Rcpp::export]]
NumericMatrix makeSigmasemiPD(NumericMatrix geno){
  Environment pkg = Environment::namespace_env("makesigma");
  Function f = pkg["makesigma"];
  return f(geno);
}

// List nearPDR(NumericMatrix geno){
//   Environment pkg = Environment::namespace_env("Matrix");
//   Function f = pkg["nearPD"];
//   return f(Named("x")=geno);
// }
// 
// //[[Rcpp::export]]
// NumericMatrix makeSigmanearPD(NumericMatrix geno){
//   Function asmatrix("as.matrix");
// 
//   mat genomat(geno.begin(), geno.nrow(), geno.ncol(), false);
//   mat raw_sigmamat = cor(genomat.t(), genomat.t());
//   NumericMatrix sigma = wrap(raw_sigmamat);
//   List res = nearPDR(sigma);
//   NumericMatrix final_res = asmatrix(res["mat"]);
//   return final_res;
// }

//[[Rcpp::export]]
DataFrame computeLRT(NumericMatrix stat, NumericMatrix geno, int subsize = 50, int causalCount = 2, double NCP = 5.7, double gamma = 0.01){

  int N = geno.ncol();
  
  NumericMatrix stat_50(subsize,1);
  NumericMatrix geno_50(subsize,N);
  
  NumericVector idx = orderC(stat);
  
  if(stat.nrow() != subsize){
    for(int i ; i <subsize; i++){
      stat_50(i,_) = stat(idx[i],_);
    }
    for(int i ; i <subsize; i++){
      geno_50(i,_) = geno(idx[i],_);
    }
  }else{
    for(int i ; i <subsize; i++){
      stat_50(i,_) = stat(i,_);
    }
    for(int i ; i <subsize; i++){
      geno_50(i,_) = geno(i,_);
    }
  }
  
  // make sigma matrix
  // NumericMatrix sigmaMatrix = makeSigmanearPD(geno_50);
  NumericMatrix sigmaMatrix = makeSigmasemiPD(geno_50);
  
  // make tot index 
  // making tot_causalIndex using R because R faster than Rcpp
  // List tot_causalIndex = combn(subsize,causalCount);
  List tot_causalIndex = combnR_marsRR(subsize,causalCount);
  
  int len_causalIndex = tot_causalIndex.size();
  
  int baseValue = 0;
  double sumLikelihood = 0;
  double allZero_likelihood = 0;
  double tmp_likelihood = 0;
  
  for(int z=0; z < len_causalIndex; z++){
    
    if( z == 0){
      int maxVal = 0;
      for(int i = 0; i < subsize; i++) {
        double s_stat_50 = stat_50(i,0);
        if (maxVal < abs(s_stat_50))
          maxVal = s_stat_50;
      }
      baseValue = maxVal * maxVal;
      allZero_likelihood = exp(-baseValue/2) * (pow(gamma, 0))*(pow(1-gamma, subsize-0));
    }
    
    NumericVector causalIndex = tot_causalIndex[z];
    
    int causalCount = causalIndex.size();
    
    mat Rcc(causalCount, causalCount, fill::zeros);
    mat Zcc(causalCount, 1, fill::zeros);
    mat mean(causalCount, 1, fill::zeros);
    mat diagC(causalCount, causalCount, fill::zeros);
    
    for (int i = 0; i < causalCount; i++){
      for(int j = 0; j < causalCount; j++) {
        Rcc(i,j) = sigmaMatrix(causalIndex[i], causalIndex[j]);
      }
      Zcc(i,0) = stat_50(causalIndex[i],0);
      diagC(i,i) = NCP;
    }
    
    while (det(Rcc) <= 0.01) {
      mat toAdd(causalCount, causalCount);
      toAdd.eye();
      Rcc = Rcc + 0.1 * toAdd;
    }
    
    tmp_likelihood = fracdmvnorm(Zcc, mean, Rcc, diagC, NCP, baseValue) * (pow(gamma, causalCount))*(pow(1-gamma, subsize-causalCount));
    sumLikelihood += tmp_likelihood;
  }
  
  double LRTscore = sumLikelihood/allZero_likelihood;//compute LRTscore
  
  NumericVector res = computePvalue(stat_50);
  res.push_front(LRTscore);
  res.names() = CharacterVector({"LRT","uni"});
  DataFrame save_res = DataFrame::create( Named("LRT") = res[0] , _["uni"] = res[1]);
  return(save_res);
  
  // res.attr("dim") = Dimension(1, 2);
  // NumericMatrix save_res = as<NumericMatrix>(res);
  // 
  // CharacterVector ch = {"LRT","uni"};
  // colnames(save_res) = ch;
  // return(save_res);
}

void set_seedC(double seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(std::floor(std::fabs(seed)));
}

// [[Rcpp::export]]
arma::mat rmvnormC(int n, arma::vec mean, arma::mat sigma, double setseed = 1) {
  int ncols = sigma.n_cols;
  
  set_seedC(setseed);
  mat Y = randn(n, ncols);
  
  return repmat(mean, 1, n).t() + Y * chol(sigma);
}


static double const log2pi = std::log(2.0 * M_PI);

void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;
  
  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

// [[Rcpp::export]]
arma::vec dmvnormC(arma::mat const &x, arma::rowvec const &mean, arma::mat const &sigma, bool const logd = false) { 
  using arma::uword;
  uword const n = x.n_rows, 
    xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())), 
    constants = -(double)xdim/2.0 * log2pi, 
    other_terms = rootisum + constants;
  
  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean);
    inplace_tri_mat_mult(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);     
  }  
  
  if (logd)
    return out;
  return exp(out);
}

// [[Rcpp::export]]
List nullmars(NumericMatrix geno, int simNum = 10000, double setseed = 1, bool fastMARS = true) {
  int n = geno.ncol();

  NumericMatrix raw_xs = scaleC(transpose(geno));
  mat xs(raw_xs.begin(), raw_xs.nrow(), raw_xs.ncol(), false);

  mat I(n,n,fill::eye);
  vec mu = zeros<vec>(n);
  
  if(fastMARS){
    mat I2 = sqrt(2) * I;
    mat Sall = rmvnormC(simNum, mu, I2, setseed);
    mat raw_Sall_new = (xs.t() / sqrt(n)) * Sall.t();
    NumericMatrix Sall_new = wrap(raw_Sall_new);
    
    rowvec row_mu = zeros<rowvec>(n);
    vec raw_w_1 = dmvnormC(Sall, row_mu, I, true);
    vec raw_w_2 = dmvnormC(Sall, row_mu, I2, true);
    vec raw_w(simNum, fill::zeros);
    for (int i = 0; i < simNum; i++){
      raw_w[i] = raw_w_1[i] / raw_w_2[i];
    }

    NumericMatrix w = wrap(raw_w);
    List nulldat = List::create(Named("w") = w, _["Sall_new"] = Sall_new);
    return(nulldat);
  }else{
    mat Sall = rmvnormC(simNum, mu, I, setseed);
    mat raw_Sall_new = (xs.t() / sqrt(n)) * Sall.t();
    NumericMatrix Sall_new = wrap(raw_Sall_new);
    NumericMatrix w (simNum, 1);
    
    List nulldat = List::create(Named("w") = w, _["Sall_new"] = Sall_new);
    return(nulldat);
  }

}
