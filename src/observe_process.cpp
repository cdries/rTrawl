#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


List observe_process(arma::vec x_grid_latent, arma::vec p_grid_latent, 
                     double T0, double TT, double observed_freq) {
  // returns time stamps and process values at a chosen observation frequency,
  // handles non-equidistant observations and does not return an equidistant grid
  //
  // arguments:
  // x_grid_latent  : vector with observed times
  // p_grid_latent  : corresponding process values
  // T0             : beginpoint of observation interval
  // TT             : endpoint of observation interval
  // observed_freq  : new observation frequency
  //
  // author: Dries Cornilly

  double Mprec = std::numeric_limits<double>::epsilon();
  
  // subset to -Inf - TT
  arma::uvec ind_TT = find(x_grid_latent < TT + Mprec);
  x_grid_latent = x_grid_latent.elem(ind_TT);
  p_grid_latent = p_grid_latent.elem(ind_TT);
  
  // create observation grid
  arma::vec x_grid_observed = T0 + arma::ceil((x_grid_latent - T0 - Mprec) / observed_freq) * observed_freq;
  int n = x_grid_observed.n_elem;
  
  arma::vec ind = arma::zeros(n);
  for (int ii = 1; ii < n; ii++) {
    if (x_grid_observed(ii) - x_grid_observed(ii - 1) > observed_freq / 2.0) ind(ii - 1) = 1;
  }
  ind(n - 1) = 1;
  
  // deal with observation at T0
  if (x_grid_latent(0) < T0 + Mprec) {
    int ii = 0;
    while (x_grid_observed(ii) < T0 + observed_freq / 2.0) {
      ind(ii) = 0;
      ii++;
    }
    ind(ii - 1) = 1;
    x_grid_observed(ii - 1) = T0;
  }
  
  // subset to observation grid
  arma::uvec ind_observed = find(ind > 0.5);
  x_grid_observed = x_grid_observed.elem(ind_observed);
  arma::vec p_grid_observed = p_grid_latent.elem(ind_observed);
  
  // return list
  List out;
  out["x_grid_observed"] = x_grid_observed;
  out["p_grid_observed"] = p_grid_observed;
  
  return out;
}
