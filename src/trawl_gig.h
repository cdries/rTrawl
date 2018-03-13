#ifndef TRAWL_GIG_H
#define TRAWL_GIG_H

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


arma::vec survival_GIG(arma::vec unif_seed, double gamma, double delta, 
                       double nu, double Tmax, double b, double observed_freq);

arma::vec leb_AtA_GIG(arma::vec h, double gamma, double delta, double nu);

arma::mat d_leb_AtA_GIG(arma::vec h, double gamma, double delta, double nu);

#endif
