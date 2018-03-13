#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


double intens_NEGBIN(double m, double theta) {
  
  double intens = -m * log(1.0 - theta);
  
  return intens;
}

arma::vec rjump_NEGBIN(int n, double theta) {
  
  arma::vec rj = runif(n, 0.0, 1.0);
  double log1theta = -log(1.0 - theta);
  
  for (int ii = 0; ii < n; ii++) {
    double uii = rj(ii);
    int k = 1;
    double cumprob = theta / log1theta;
    while (uii > cumprob) {
      k++;
      cumprob += pow(theta, k) / (k * log1theta);
    }
    rj(ii) = k;
  }
  
  return rj;
}

double cum_NEGBIN(int ord, double m, double theta) {
  
  double cum = 0.0;
  if (ord == 1) cum = m * theta / (1 - theta);
  if (ord == 2) cum = m * theta / ((1 - theta) * (1 - theta));
  
  return cum;
}
