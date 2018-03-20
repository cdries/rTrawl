#ifndef CUM_H
#define CUM_H

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


double cum_sample(int ord, arma::vec x_grid, arma::vec p_grid, double TT);

#endif
