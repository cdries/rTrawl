#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

double intens_SKELLAM(double nu_p, double nu_m) {
  // Poisson intensity for the Skellam levy process,
  // see Barndorff-Nielsen, Lunde, Shephard and Veraart (2014)
  //
  // arguments:
  // nu_p : parameter of the Skellam process, positive part
  // nu_m : parameter of the Skellam process, negative part
  //
  // author: Dries Cornilly
  
  double intens = nu_p + nu_m;
  
  return intens;
}

arma::vec rjump_SKELLAM(int n, double nu_p, double nu_m) {
  // sample jump sizes for the Skellam levy process,
  // see Barndorff-Nielsen, Lunde, Shephard and Veraart (2014)
  //
  // arguments:
  // n    : sample size
  // nu_p : parameter of the Skellam process, positive part
  // nu_m : parameter of the Skellam process, negative part
  //
  // author: Dries Cornilly
  
  arma::vec rj = 2.0 * (rbinom(n, 1, nu_p / (nu_p + nu_m)) - 0.5);
  
  return rj;
}

arma::vec fit_SKELLAM(double k1_sample, double k2_sample) {
  // fit Skellam distribution based on its first two cumulants
  // see Barndorff-Nielsen, Lunde, Shephard and Veraart (2014)
  //
  // arguments:
  // k1_sample  : first cumulant of the Skellam distribution
  // k2_sample  : second cumulant of the Skellam distribution
  //
  // author: Dries Cornilly
  
  arma::vec levy_par = arma::ones(2);
  levy_par(0) = 0.5 * (k1_sample + k2_sample);
  levy_par(1) = 0.5 * (k2_sample - k1_sample);
  
  return levy_par;
}
