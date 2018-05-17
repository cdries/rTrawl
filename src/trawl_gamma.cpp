#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


arma::vec survival_GAMMA(arma::vec unif_seed, double alpha, double H, double Tmax, double b) {
  // computes the time an observation is covered by the gamma trawl
  //
  // arguments:
  // unif_seed  : standard uniform randomly generated numbers
  // alpha      : trawl parameter
  // H          : trawl parameter
  // Tmax       : maximum observation window, times larger than this are not useful
  // b          : b parameters, as in Shephard and Yang (2017)
  //
  // author: Dries Cornilly
  
  arma::vec st = Tmax * arma::ones(unif_seed.n_elem);
  arma::uvec ind_trawl = find(unif_seed > b);
  st.elem(ind_trawl) = alpha * (arma::pow((unif_seed.elem(ind_trawl) - b) / (1.0 - b), -1.0 / H) - 1.0);
  
  return st;
}

arma::vec leb_AtA_GAMMA(arma::vec h, double alpha, double H) {
  // computes lebesgue measure of A_0 cap A_h for the gamma trawl
  //
  // arguments:
  // h          : determines A_h, if 0, then leb(A) is computed
  // alpha      : trawl parameter
  // H          : trawl parameter
  //
  // author: Dries Cornilly
  
  arma::vec leb = alpha * arma::pow(1.0 + h / alpha, 1.0 - H) / (H - 1.0);
  
  return leb;
}

arma::mat d_leb_AtA_GAMMA(arma::vec h, double alpha, double H) {
  // derivative of the lebesgue measure of A_0 cap A_h with respect to the 
  // trawl parameters
  //
  // arguments:
  // h          : determines A_h, if 0, then leb(A) is computed
  // alpha      : trawl parameter
  // H          : trawl parameter
  //
  // author: Dries Cornilly
  
  arma::vec ha1 = 1.0 + h / alpha;
  arma::vec ha11H = arma::pow(ha1, 1.0 - H);
  
  arma::mat d_leb = arma::zeros(h.n_elem, 2);
  d_leb.col(0) = (ha11H - (1.0 - H) * h % arma::pow(ha1, -H) / alpha) / (H - 1.0);
  d_leb.col(1) = -alpha * ((H - 1.0) * log(ha1) % ha11H + ha11H) / ((H - 1.0) * (H - 1.0));
  
  return d_leb;
}

arma::vec trawl_GAMMA(arma::vec h, double alpha, double H) {
  // gamma trawl function
  //
  // arguments:
  // h          : argument at which to compute the function
  // alpha      : trawl parameter
  // H          : trawl parameter
  //
  // author: Dries Cornilly
  
  arma::vec val = arma::zeros(h.n_elem);
  arma::uvec ind = find(h < std::numeric_limits<double>::epsilon());
  val.elem(ind) = arma::pow(1.0 - h.elem(ind) / alpha, -H);
  
  return val;
}
