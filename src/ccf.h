#ifndef CCF_H
#define CCF_H

#include "RcppArmadillo.h"
#include "trawl_wrap.h"
#include "observe_process.h"
#include "cum.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


double ccf_helper(arma::vec x1_grid, arma::vec p1_grid, arma::vec x2_grid, 
                  arma::vec p2_grid, double TT, double h);

#endif
