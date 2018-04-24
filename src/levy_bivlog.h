#ifndef LEVY_BIVLOG_H
#define LEVY_BIVLOG_H

#include "RcppArmadillo.h"
#include "levy_negbin.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


arma::mat rjump_BIVLOG(int n, double p1, double p2);

double intens_BIVLOG(double m, double theta1, double theta2);

#endif
