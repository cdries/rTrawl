#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


arma::vec survival_EXP(arma::vec unif_seed, double lambda, double Tmax, double b) {
  
  arma::vec st = Tmax * arma::ones(unif_seed.n_elem);
  arma::uvec ind_trawl = find(unif_seed > b);
  st.elem(ind_trawl) = -1.0 * log((unif_seed.elem(ind_trawl) - b) / (1.0 - b)) / lambda;
  
  return st;
}

arma::vec leb_AtA_EXP(arma::vec h, double lambda) {

  arma::vec leb = exp(-lambda * h) / lambda;

  return leb;
}
