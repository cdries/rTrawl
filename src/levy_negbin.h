#ifndef LEVY_NEGBIN_H
#define LEVY_NEGBIN_H

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


double intens_NEGBIN(double m, double theta);

arma::vec rjump_NEGBIN(int n, double theta);

double cum_NEGBIN(int ord, double m, double theta);

#endif
