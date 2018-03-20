#ifndef LEVY_WRAP_H
#define LEVY_WRAP_H

#include "RcppArmadillo.h"
#include "levy_poisson.h"
#include "levy_skellam.h"
#include "levy_negbin.h"
#include "levy_dnegbin.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

double levy_intens(std::string levy_seed, arma::vec levy_par);

arma::vec levy_rjump(int n, std::string levy_seed, arma::vec levy_par);

arma::vec levy_cum_fit(std::string levy_seed, double k1_sample, double k2_sample);

arma::mat levy_varcovar(std::string levy_seed, arma::mat levy_par, arma::mat design_matrix);

#endif
