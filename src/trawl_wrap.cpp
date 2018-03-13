#include "RcppArmadillo.h"
#include "trawl_exp.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


arma::vec trawl_times(arma::vec unif_seed, std::string trawl, arma::vec trawl_par,
                      double observed_freq, double Tmax, double b) {
  
  arma::vec surv_times(unif_seed.n_elem);
  if (trawl == "exp") {
    surv_times = survival_EXP(unif_seed, trawl_par(0), Tmax, b);
  // } else if (trawl == "gamma") {
  //   surv_times = survival_GAMMA(unif_seed, trawl_par(0), trawl_par(1), Tmax, b);
  // } else if (trawl == "invGauss") {
  //   surv_times = survival_INVGAUSS(unif_seed, trawl_par(0), trawl_par(1), Tmax, b, observed_freq);
  // } else if (trawl == "gig") {
  //   surv_times = survival_GIG(unif_seed, trawl_par(0), trawl_par(1), trawl_par(2), Tmax, b, observed_freq);
  } else {
    stop("provide a valid trawl");
  }
  
  return surv_times;
}
