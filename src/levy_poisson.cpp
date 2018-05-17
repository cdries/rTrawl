#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

double intens_POISSON(double nu) {
  // Poisson intensity for the Poisson levy process,
  // see Barndorff-Nielsen, Lunde, Shephard and Veraart (2014)
  //
  // arguments:
  // nu   : parameter of the Poisson distribution
  //
  // author: Dries Cornilly
  
  double intens = nu;
  
  return intens;
}

arma::vec rjump_POISSON(int n) {
  // sample jump sizes for the Poisson levy process,
  // see Barndorff-Nielsen, Lunde, Shephard and Veraart (2014)
  //
  // arguments:
  // n    : sample size
  //
  // author: Dries Cornilly
  
  arma::vec rj = arma::ones(n);
  
  return rj;
}

arma::vec fit_POISSON(double k1_sample) {
  // fit Poisson distribution based on its first moment
  // see Barndorff-Nielsen, Lunde, Shephard and Veraart (2014)
  //
  // arguments:
  // k1_sample  : first cumulant of the Poisson distribution
  //
  // author: Dries Cornilly
  
  arma::vec levy_par = arma::ones(1) * k1_sample;
  
  return levy_par;
}
