#include "RcppArmadillo.h"
#include "levy_negbin.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


arma::mat rjump_BIVLOG(int n, double p1, double p2) {
  
  arma::mat rj = arma::zeros(n, 2);
  rj.col(0) = rjump_NEGBIN(n, p1 / (1.0 - p2));
  
  double delta1 = log(1.0 - p2) / log(1.0 - p1 - p2);
  rj.col(0) *= rbinom(n, 1, 1.0 - delta1)(0);
  
  for (int ii = 0; ii < n; ii++) {
    if (rj(ii, 0) < 0.5) {
      rj(ii, 1) = rjump_NEGBIN(n, p2)(0);
    } else {
      rj(ii, 1) = rnbinom(1, rj(ii, 0), p2)(0);
    }
  }
  
  return rj;
}

double intens_BIVLOG(double m, double theta1, double theta2) {
  
  double alpha1 = 1.0 / (1.0 + theta1);
  double alpha2 = 1.0 / (1.0 + theta2);
  double intens = m * log(1.0 + alpha1 + alpha2);
  
  return intens;
}
