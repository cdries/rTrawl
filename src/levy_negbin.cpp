#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


double intens_NEGBIN(double m, double theta) {
  // Poisson intensity for the negative binomial levy process,
  // see Barndorff-Nielsen, Lunde, Shephard and Veraart (2014)
  //
  // arguments:
  // m        : first parameter in the negative binomial
  // theta    : second parameter in the negative binomial
  //
  // author: Dries Cornilly
  
  double intens = -m * log(1.0 - theta);
  
  return intens;
}

arma::vec rjump_NEGBIN(int n, double theta) {
  // sample jump sizes for the negative binomial levy process,
  // see Barndorff-Nielsen, Lunde, Shephard and Veraart (2014)
  //
  // arguments:
  // n        : sample size
  // theta    : second parameter in the negative binomial
  //
  // author: Dries Cornilly
  
  arma::vec rj = runif(n, 0.0, 1.0);
  double log1theta = -log(1.0 - theta);
  
  for (int ii = 0; ii < n; ii++) {
    double uii = rj(ii);
    int k = 1;
    double cumprob = theta / log1theta;
    while (uii > cumprob) {
      k++;
      cumprob += pow(theta, k) / (k * log1theta);
    }
    rj(ii) = k;
  }
  
  return rj;
}

double cum_NEGBIN(int ord, double m, double theta) {
  // compute cumulants of the negative binomial distribution,
  // see Barndorff-Nielsen, Lunde, Shephard and Veraart (2014)
  //
  // arguments:
  // ord      : order of the cumulant - 1 or 2
  // m        : first parameter in the negative binomial
  // theta    : second parameter in the negative binomial
  //
  // author: Dries Cornilly
  
  double cum = 0.0;
  if (ord == 1) cum = m * theta / (1 - theta);
  if (ord == 2) cum = m * theta / ((1 - theta) * (1 - theta));
  
  return cum;
}

arma::vec fit_NEGBIN(double k1_sample, double k2_sample) {
  // fit Negative Binomial distribution based on its first two cumulants
  // see Barndorff-Nielsen, Lunde, Shephard and Veraart (2014)
  //
  // arguments:
  // k1_sample  : first cumulant of the Skellam distribution
  // k2_sample  : second cumulant of the Skellam distribution
  //
  // author: Dries Cornilly
  
  arma::vec levy_par = arma::ones(2);
  levy_par(1) = 1.0 - k1_sample / k2_sample;
  if (levy_par(1) < 0) levy_par(1) = 0.001;
  levy_par(0) = (1.0 - levy_par(1)) * k1_sample / levy_par(1);
  
  return levy_par;
}
