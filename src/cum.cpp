#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export()]]
double cum_sample(int ord, arma::vec x_grid, arma::vec p_grid, double T0, double TT) {
  
  // mean and centered observations
  int n = x_grid.n_elem;
  arma::vec diff_x = arma::diff(x_grid);
  double k1L = (arma::sum(diff_x % p_grid.head(n - 1)) +
                (TT - x_grid(n - 1)) * p_grid(n - 1)) / (TT - p_grid(0));
  
  double cum = 0.0;
  if (ord == 1) {
    cum = k1L;
    
  } else if (ord == 2) {
    p_grid -= k1L;
    
    // variance
    cum = (arma::sum(diff_x % p_grid.head(n - 1) % p_grid.head(n - 1)) +
                  (TT - x_grid(n - 1)) * p_grid(n - 1) * p_grid(n - 1)) / (TT - p_grid(0));
  }
  
  return cum;
}
