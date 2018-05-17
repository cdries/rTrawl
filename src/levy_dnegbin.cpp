#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


double intens_DNEGBIN(double m_p, double theta_p, double m_m, double theta_m) {
  // Poisson intensity for the difference of two negative binomial levy processes,
  // see Barndorff-Nielsen, Lunde, Shephard and Veraart (2014)
  //
  // arguments:
  // m_p      : first parameter in the negative binomial, positive component
  // theta_p  : second parameter in the negative binomial, positive component
  // m_m      : first parameter in the negative binomial, negative component
  // theta_m  : second parameter in the negative binomial, negative component
  //
  // author: Dries Cornilly
  
  double intens = -m_p * log(1.0 - theta_p) - m_m * log(1.0 - theta_m);
  
  return intens;
}

arma::vec rjump_DNEGBIN(int n, double m_p, double theta_p, double m_m, double theta_m) {
  // sample jump sizes for the difference of two negative binomial levy processes,
  // see Barndorff-Nielsen, Lunde, Shephard and Veraart (2014)
  //
  // arguments:
  // n        : number of observations
  // m_p      : first parameter in the negative binomial, positive component
  // theta_p  : second parameter in the negative binomial, positive component
  // m_m      : first parameter in the negative binomial, negative component
  // theta_m  : second parameter in the negative binomial, negative component
  //
  // author: Dries Cornilly
  
  arma::vec rj = runif(n, 0.0, 1.0);
  double log1theta_p = -log(1.0 - theta_p);
  double log1theta_m = -log(1.0 - theta_m);
  double prob_pos = m_p * log1theta_p / (m_p * log1theta_p + m_m * log1theta_m);
  arma::vec posneg = runif(n, 0.0, 1.0);
  
  for (int ii = 0; ii < n; ii++) {
    
    double uii = rj(ii);
    int k = 1;
    
    if (posneg(ii) < prob_pos) {
      double cumprob = theta_p / log1theta_p;
      while (uii > cumprob) {
        k++;
        cumprob += pow(theta_p, k) / (k * log1theta_p);
      }
      rj(ii) = k;
    } else {
      double cumprob = theta_m / log1theta_m;
      while (uii > cumprob) {
        k++;
        cumprob += pow(theta_m, k) / (k * log1theta_m);
      }
      rj(ii) = -k;
    }
  }
  
  return rj;
}

double cum_DNEGBIN(int ord, double m_p, double theta_p, double m_m, double theta_m) {
  // compute cumulants for the distribution of the difference of two negative binomials,
  // see Barndorff-Nielsen, Lunde, Shephard and Veraart (2014)
  //
  // arguments:
  // ord      : order of the cumulant - 1 or 2
  // m_p      : first parameter in the negative binomial, positive component
  // theta_p  : second parameter in the negative binomial, positive component
  // m_m      : first parameter in the negative binomial, negative component
  // theta_m  : second parameter in the negative binomial, negative component
  //
  // author: Dries Cornilly
  
  double cum = 0.0;
  if (ord == 1) cum = m_p * theta_p / (1 - theta_p) - m_m * theta_m / (1 - theta_m);
  if (ord == 2) cum = m_p * theta_p / ((1 - theta_p) * (1 - theta_p)) +
    m_m * theta_m / ((1 - theta_m) * (1 - theta_m));
  
  return cum;
}
