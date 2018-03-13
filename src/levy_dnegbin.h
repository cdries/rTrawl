#ifndef LEVY_DNEGBIN_H
#define LEVY_DNEGBIN_H

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


double intens_DNEGBIN(double m_p, double theta_p, double m_m, double theta_m);

arma::vec rjump_DNEGBIN(int n, double m_p, double theta_p, double m_m, double theta_m);

double cum_DNEGBIN(int ord, double m_p, double theta_p, double m_m, double theta_m);

#endif
