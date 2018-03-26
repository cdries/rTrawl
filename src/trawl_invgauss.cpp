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

arma::vec leb_AtA_INVGAUSS(arma::vec h, double gamma, double delta) {
  
  arma::vec leb = gamma / delta *
    exp(delta * gamma * (1.0 - arma::sqrt(1.0 + 2.0 * h / (gamma * gamma))));
  
  return leb;
}

arma::mat d_leb_AtA_INVGAUSS(arma::vec h, double gamma, double delta) {
  
  double gamma2 = gamma * gamma;
  arma::vec sqrtdg = arma::sqrt(1.0 + 2.0 * h / gamma2);
  arma::vec sqrtdg_inv = arma::pow(sqrtdg, -1.0);
  arma::vec dg1s = delta * gamma * (1.0 - sqrtdg);
  
  arma::mat d_leb = arma::zeros(h.n_elem, 2);
  d_leb.col(0) = exp(dg1s) / delta % 
    (1.0 + gamma * delta * (1.0 - sqrtdg + 2.0 * h % sqrtdg_inv / gamma2));
  d_leb.col(1) = gamma * exp(dg1s) % (dg1s - 1.0) / (delta * delta);
  
  return d_leb;
}

arma::vec trawl_INVGAUSS(arma::vec h, double gamma, double delta) {
  
  arma::vec val = arma::zeros(h.n_elem);
  arma::uvec ind = find(h < std::numeric_limits<double>::epsilon());
  
  arma::vec z = arma::sqrt(1.0 - 2.0 * h.elem(ind) / (gamma * gamma));
  val.elem(ind) = arma::pow(z, -1.0) * exp(delta * gamma * (1.0 - z));
  
  return val;
}
