#ifndef OBSERVE_PROCESS_H
#define OBSERVE_PROCESS_H

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

List observe_process(arma::vec x_grid_latent, arma::vec p_grid_latent, 
                     double T0, double TT, double observed_freq);

#endif
