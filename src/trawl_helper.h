#ifndef TRAWL_HELPER_H
#define TRAWL_HELPER_H

#include "RcppArmadillo.h"
#include "trawl_wrap.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


double leb_EXP_EXP(double h, double lambda1, double lambda2);

double leb_GEN_GEN(double  h, std::string trawl1, arma::vec trawl_par1, 
                   std::string trawl2, arma::vec trawl_par2);

#endif
