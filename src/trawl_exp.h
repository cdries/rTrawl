#ifndef TRAWL_EXP_H
#define TRAWL_EXP_H

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


arma::vec survival_EXP(arma::vec unif_seed, double lambda, double Tmax, double b);

arma::vec leb_AtA_EXP(arma::vec h, double lambda);

arma::mat d_leb_AtA_EXP(arma::vec h, double lambda);

arma::vec trawl_EXP(arma::vec h, double lambda);

#endif
