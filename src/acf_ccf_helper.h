#ifndef ACF_CCF_HELPER_H
#define ACF_CCF_HELPER_H

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


double ccf_helper(arma::vec x1_grid, arma::vec p1_grid, arma::vec x2_grid, 
                  arma::vec p2_grid, double TT, double h);

double ccf_dp_helper(double h, arma::vec x_grid1, arma::vec diff_p1, arma::vec x_grid2, 
                     arma::vec diff_p2, double TT, int lg);

#endif
