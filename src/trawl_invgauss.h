#ifndef TRAWL_INVGAUSS_H
#define TRAWL_INVGAUSS_H

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


arma::vec survival_INVGAUSS(arma::vec unif_seed, double gamma, double delta, 
                            double Tmax, double b, double observed_freq);

arma::vec leb_AtA_INVGAUSS(arma::vec h, double gamma, double delta);

arma::mat d_leb_AtA_INVGAUSS(arma::vec h, double gamma, double delta);

arma::vec trawl_INVGAUSS(arma::vec h, double gamma, double delta);

#endif
