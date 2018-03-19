#ifndef LEVY_SKELLAM_H
#define LEVY_SKELLAM_H

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

double intens_SKELLAM(double nu_p, double nu_m);

arma::vec rjump_SKELLAM(int n, double nu_p, double nu_m);

arma::vec fit_SKELLAM(double k1_sample, double k2_sample);

#endif
