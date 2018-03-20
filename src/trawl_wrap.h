#ifndef TRAWL_WRAP_H
#define TRAWL_WRAP_H

#include "RcppArmadillo.h"
#include "trawl_exp.h"
#include "trawl_gamma.h"
#include "trawl_invgauss.h"
#include "trawl_gig.h"
#include "trawl_helper.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


arma::vec trawl_times(arma::vec unif_seed, std::string trawl, arma::vec trawl_par,
                      double observed_freq, double Tmax, double b);

int number_parameters_trawl(std::string trawl);

List trawl_bounds(std::string trawl);

arma::vec trawl_x0(std::string trawl);

arma::vec leb_AtA(arma::vec h, std::string trawl, arma::vec trawl_par);

arma::mat d_leb_AtA(arma::vec h, std::string trawl, arma::vec trawl_par);

arma::vec leb_autocorrelator(arma::vec h, std::string trawl1, arma::vec trawl_par1,
                             std::string trawl2, arma::vec trawl_par2);

#endif
