#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export()]]
arma::vec survival_GIG(arma::vec unif_seed, double gamma, double delta, 
                       double nu, double Tmax, double b, double observed_freq) {
  
  arma::vec st = Tmax * arma::ones(unif_seed.n_elem);
  arma::uvec ind_trawl = find(unif_seed > b);
  arma::vec u_seed = (unif_seed.elem(ind_trawl) - b) / (1.0 - b);
  
  double gd = gamma * delta;
  double Kngd = Rf_bessel_k(gd, nu, 1.0);
  double gamma2 = gamma * gamma;
  double eps = observed_freq / 10;
  
  for (int ii = 0; ii < ind_trawl.n_elem; ii++) {
    
    double uii = u_seed(ii);
    double lower = 0.0;
    double upper = Tmax;
    
    if (pow(1.0 + 2.0 * upper / gamma2, -0.5 * nu) *
        Rf_bessel_k(gd * sqrt(1.0 + 2.0 * upper / gamma2), nu, 1.0) /
          Kngd > uii) {
      st(ind_trawl(ii)) = Tmax;
    } else {
      while (upper - lower > eps) {
        double ulm = 0.5 * (upper + lower);
        if (pow(1.0 + 2.0 * ulm / gamma2, -0.5 * nu) *
            Rf_bessel_k(gd * sqrt(1.0 + 2.0 * ulm / gamma2), nu, 1.0) /
              Kngd < uii) {
          upper = ulm;
        } else {
          lower = ulm;
        }
      }
      st(ind_trawl(ii)) = lower;
    }
  }
  
  return st;
}
