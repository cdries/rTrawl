#ifndef TRAWL_GAMMA_H
#define TRAWL_GAMMA_H

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


arma::vec survival_GAMMA(arma::vec unif_seed, double alpha, double H, double Tmax, double b);

#endif
