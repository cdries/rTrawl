#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


arma::vec survival_EXP(arma::vec unif_seed, double lambda, double Tmax, double b) {
  // computes the time an observation is covered by the exponential trawl
  //
  // arguments:
  // unif_seed  : standard uniform randomly generated numbers
  // lambda     : trawl parameter
  // Tmax       : maximum observation window, times larger than this are not useful
  // b          : b parameters, as in Shephard and Yang (2017)
  //
  // author: Dries Cornilly
  
  arma::vec st = Tmax * arma::ones(unif_seed.n_elem);
  arma::uvec ind_trawl = find(unif_seed > b);
  st.elem(ind_trawl) = -1.0 * log((unif_seed.elem(ind_trawl) - b) / (1.0 - b)) / lambda;
  
  return st;
}

arma::vec leb_AtA_EXP(arma::vec h, double lambda) {
  // computes lebesgue measure of A_0 cap A_h for the exponential trawl
  //
  // arguments:
  // h          : determines A_h, if 0, then leb(A) is computed
  // lambda     : trawl parameter
  //
  // author: Dries Cornilly

  arma::vec leb = exp(-lambda * h) / lambda;

  return leb;
}

arma::mat d_leb_AtA_EXP(arma::vec h, double lambda) {
  // derivative of the lebesgue measure of A_0 cap A_h with respect to the 
  // trawl parameters
  //
  // arguments:
  // h          : determines A_h, if 0, then leb(A) is computed
  // lambda     : trawl parameter
  //
  // author: Dries Cornilly
  
  arma::mat d_leb = -exp(-lambda * h) % (lambda * h + 1.0) / (lambda * lambda);
  
  return d_leb;
}

arma::vec trawl_EXP(arma::vec h, double lambda) {
  // exponential trawl function
  //
  // arguments:
  // h          : argument at which to compute the function
  // lambda     : trawl parameter
  //
  // author: Dries Cornilly
  
  arma::vec val = arma::zeros(h.n_elem);
  arma::uvec ind = find(h < std::numeric_limits<double>::epsilon());
  val.elem(ind) = exp(lambda * h.elem(ind));
  
  return val;
}
