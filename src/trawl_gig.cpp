#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


arma::vec survival_GIG(arma::vec unif_seed, double gamma, double delta, 
                       double nu, double Tmax, double b, double observed_freq) {
  // computes the time an observation is covered by the gig trawl
  //
  // arguments:
  // unif_seed  : standard uniform randomly generated numbers
  // gamma      : trawl parameter
  // delta      : trawl parameter
  // nu         : trawl parameter
  // Tmax       : maximum observation window, times larger than this are not useful
  // b          : b parameters, as in Shephard and Yang (2017)
  // observed_freq : determines how accurate survival times should be approximated
  //
  // author: Dries Cornilly
  
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
      
      // survival time is longer than maximum observation time
      st(ind_trawl(ii)) = Tmax;
      
    } else {
      
      // divide and search algorithm untill the interval length containing
      // the survival time is within the required accuracy
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

arma::vec leb_AtA_GIG(arma::vec h, double gamma, double delta, double nu) {
  // computes lebesgue measure of A_0 cap A_h for the gig trawl
  //
  // arguments:
  // h          : determines A_h, if 0, then leb(A) is computed
  // gamma      : trawl parameter
  // delta      : trawl parameter
  // nu         : trawl parameter
  //
  // author: Dries Cornilly
  
  arma::vec z = arma::sqrt(1.0 + 2.0 * h / (gamma * gamma));
  double besseldg = Rf_bessel_k(delta * gamma, nu, 1.0);
  int n_h = h.n_elem;
  
  arma::vec leb = arma::zeros(n_h);
  for (int ii = 0; ii < n_h; ii++) {
    leb(ii) = gamma / delta * pow(z(ii), 1.0 - nu) *
      Rf_bessel_k(delta * gamma * z(ii), nu - 1.0, 1.0) / besseldg;
  }

  return leb;
}

arma::mat d_leb_AtA_GIG(arma::vec h, double gamma, double delta, double nu) {
  // derivative of the lebesgue measure of A_0 cap A_h with respect to the 
  // trawl parameters
  //
  // arguments:
  // h          : determines A_h, if 0, then leb(A) is computed
  // gamma      : trawl parameter
  // delta      : trawl parameter
  // nu         : trawl parameter
  //
  // author: Dries Cornilly
  
  double eps = 1e-7;
  arma::vec leb_base = leb_AtA_GIG(h, gamma, delta, nu);
  
  arma::mat d_leb = arma::zeros(h.n_elem, 3);
  d_leb.col(0) = (leb_AtA_GIG(h, gamma + eps, delta, nu) - leb_base) / eps;
  d_leb.col(1) = (leb_AtA_GIG(h, gamma, delta + eps, nu) - leb_base) / eps;
  d_leb.col(2) = (leb_AtA_GIG(h, gamma, delta, nu + eps) - leb_base) / eps;
  
  return d_leb;
}

arma::vec trawl_GIG(arma::vec h, double gamma, double delta, double nu) {
  // gig trawl function
  //
  // arguments:
  // h          : argument at which to compute the function
  // gamma      : trawl parameter
  // delta      : trawl parameter
  // nu         : trawl parameter
  //
  // author: Dries Cornilly
  
  int n_h = h.n_elem;
  arma::vec z = arma::sqrt(1.0 - 2.0 * h / (gamma * gamma));
  arma::vec val = arma::zeros(n_h);
  for (int ii = 0; ii < n_h; ii++) {
    if (h(ii) < std::numeric_limits<double>::epsilon()) {
      val(ii) = Rf_bessel_k(gamma * delta * z(ii), nu, 1.0);
    }
  }
  val %= arma::pow(z, -nu) / Rf_bessel_k(gamma * delta, nu, 1.0);

  return val;
}
