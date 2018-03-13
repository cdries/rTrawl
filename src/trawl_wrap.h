#ifndef TRAWL_WRAP_H
#define TRAWL_WRAP_H

#include "RcppArmadillo.h"
#include "trawl_exp.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

arma::vec trawl_times(arma::vec unif_seed, std::string trawl, arma::vec trawl_par,
                      double observed_freq, double Tmax, double b);

#endif
