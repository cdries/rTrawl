#ifndef CUM_H
#define CUM_H

#include "RcppArmadillo.h"
#include "trawl_wrap.h"
#include "acf_ccf_helper.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


double cum_sample(int ord, arma::vec x_grid, arma::vec p_grid, double TT);

double cum_dp_sample(int ord, double h, arma::vec x_grid, arma::vec diff_p, double TT);

List levy_cum_mv2fit(double T0, double TT, List x_grid, List p_grid, List trawl, List trawl_par, int p);

#endif
