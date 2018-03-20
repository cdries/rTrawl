#ifndef CUM_H
#define CUM_H

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


double cum_sample(int ord, arma::vec x_grid, arma::vec p_grid, double TT);

double cum_dp_sample(int ord, double h, arma::vec x_grid, arma::vec diff_p, double TT);

#endif
