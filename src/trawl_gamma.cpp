#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


arma::vec survival_GAMMA(arma::vec unif_seed, double alpha, double H, double Tmax, double b) {
  
  arma::vec st = Tmax * arma::ones(unif_seed.n_elem);
  arma::uvec ind_trawl = find(unif_seed > b);
  st.elem(ind_trawl) = alpha * (arma::pow((unif_seed.elem(ind_trawl) - b) / (1.0 - b), -1.0 / H) - 1.0);
  
  return st;
}

arma::vec leb_AtA_GAMMA(arma::vec h, double alpha, double H) {
  
  arma::vec leb = alpha * arma::pow(1.0 + h / alpha, 1.0 - H) / (H - 1.0);
  
  return leb;
}

arma::mat d_leb_AtA_GAMMA(arma::vec h, double alpha, double H) {
  
  arma::vec ha1 = 1.0 + h / alpha;
  arma::vec ha11H = arma::pow(ha1, 1.0 - H);
  
  arma::mat d_leb = arma::zeros(h.n_elem, 2);
  d_leb.col(0) = (ha11H - (1.0 - H) * h % arma::pow(ha1, -H) / alpha) / (H - 1.0);
  d_leb.col(1) = -alpha * ((H - 1.0) * log(ha1) % ha11H + ha11H) / ((H - 1.0) * (H - 1.0));
  
  return d_leb;
}

arma::vec trawl_GAMMA(arma::vec h, double alpha, double H) {
  
  arma::vec val = arma::pow(1.0 - h / alpha, -H);
  
  return val;
}
