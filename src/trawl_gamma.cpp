#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


arma::vec survival_GAMMA(arma::vec unif_seed, double alpha, double H, double Tmax, double b) {
  
  arma::vec st = Tmax * arma::ones(unif_seed.n_elem);
  arma::uvec ind_trawl = find(unif_seed > b);
  st.elem(ind_trawl) = alpha * (arma::pow((unif_seed.elem(ind_trawl) - b) / (1.0 - b), -1.0 / H) - 1.0);
  
  return st;
}
