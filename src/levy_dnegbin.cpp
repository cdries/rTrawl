#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


double intens_DNEGBIN(double m_p, double theta_p, double m_m, double theta_m) {
  
  double intens = -m_p * log(1.0 - theta_p) - m_m * log(1.0 - theta_m);
  
  return intens;
}

arma::vec rjump_DNEGBIN(int n, double m_p, double theta_p, double m_m, double theta_m) {
  
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
  
  double cum = 0.0;
  if (ord == 1) cum = m_p * theta_p / (1 - theta_p) - m_m * theta_m / (1 - theta_m);
  if (ord == 2) cum = m_p * theta_p / ((1 - theta_p) * (1 - theta_p)) +
    m_m * theta_m / ((1 - theta_m) * (1 - theta_m));
  
  return cum;
}
