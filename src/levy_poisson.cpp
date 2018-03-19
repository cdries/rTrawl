#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

double intens_POISSON(double nu) {
  
  double intens = nu;
  
  return intens;
}

arma::vec rjump_POISSON(int n) {
  
  arma::vec rj = arma::ones(n);
  
  return rj;
}

arma::vec fit_POISSON(double k1_sample) {
  
  arma::vec levy_par = arma::ones(1) * k1_sample;
  
  return levy_par;
}
