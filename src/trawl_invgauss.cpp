#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


arma::vec survival_INVGAUSS(arma::vec unif_seed, double gamma, double delta, 
                            double Tmax, double b, double observed_freq) {
  
  arma::vec st = Tmax * arma::ones(unif_seed.n_elem);
  arma::uvec ind_trawl = find(unif_seed > b);
  arma::vec u_seed = (unif_seed.elem(ind_trawl) - b) / (1.0 - b);
  
  double gd = gamma * delta;
  double gamma2 = gamma * gamma;
  double eps = observed_freq / 10;
  
  for (int ii = 0; ii < ind_trawl.n_elem; ii++) {
    
    double uii = u_seed(ii);
    double lower = 0.0;
    double upper = Tmax;
    
    if (pow(1.0 + 2.0 * upper / gamma2, -0.5) *
        exp(gd * (1.0 - sqrt(1.0 + 2.0 * upper / gamma2))) > uii) {
      st(ind_trawl(ii)) = Tmax;
    } else {
      while (upper - lower > eps) {
        double ulm = 0.5 * (upper + lower);
        if (pow(1.0 + 2.0 * ulm / gamma2, -0.5) *
            exp(gd * (1.0 - sqrt(1.0 + 2.0 * ulm / gamma2))) < uii) {
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
