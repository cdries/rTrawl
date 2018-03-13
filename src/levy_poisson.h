#ifndef LEVY_POISSON_H
#define LEVY_POISSON_H

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


double intens_POISSON(double nu);

arma::vec rjump_POISSON(int n);

#endif
