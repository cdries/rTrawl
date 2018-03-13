#ifndef TRAWL_EXP_H
#define TRAWL_EXP_H

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

arma::vec survival_EXP(arma::vec unif_seed, double trawl_par, double Tmax, double b);

#endif
