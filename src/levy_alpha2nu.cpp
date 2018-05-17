#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export()]]
arma::mat levy_alpha2nu(arma::mat levy_alpha, double b, double beta_0) {
  // computes nonparametric levy measure based on instantaneous jump probabilities
  // and jump rate, as in Shephard and Yang (2017)
  // 
  // arguments:
  // levy_alpha : instantaneous jump probabilities
  // b          : b parameter as in Shephard and Yang (2017)
  // beta_0     : observed jump rate
  //
  // author: Dries Cornilly
  
  // get unique integer basis
  arma::uvec indpos = arma::find(levy_alpha.col(0) % levy_alpha.col(0) > 0.5);
  arma::vec basis = levy_alpha.col(0);
  basis = unique(basis(indpos) % basis(indpos));
  basis = arma::sort(join_cols(-arma::sqrt(basis), arma::sqrt(basis)));
  
  // make basis symmetric - include zeros where necessary
  arma::vec alpha_clean = arma::zeros(basis.n_elem);
  for (int ii = 0; ii < levy_alpha.n_rows; ii++) {
    arma::uvec iipos = arma::find(basis == levy_alpha(ii, 0));
    alpha_clean(iipos(0)) = levy_alpha(ii, 1);
  }
  
  int k = alpha_clean.n_elem;
  arma::vec nu = arma::zeros(k);
  
  if (b < 0) {
    // not identifiable, assume symmetry
    for (int ii = 0; ii < k / 2; ii++) {
      double alpha_y = alpha_clean(ii);
      double alpha_miny = alpha_clean(k - ii - 1);
      nu(ii) = beta_0 * (alpha_y + alpha_miny) / 4;
    }
    
  } else {
    for (int ii = 0; ii < k / 2; ii++) {
      
      // compute levy measure
      nu(ii) = beta_0 * (alpha_clean(ii) - (1 - b) * alpha_clean(k - ii - 1)) / ((2.0 - b) * b);
      nu(k - ii - 1) = beta_0 * (alpha_clean(k - ii - 1) - (1 - b) * alpha_clean(ii)) / ((2.0 - b) * b);
      
      // make positive if required
      if (nu(ii) < 0.0) {
        nu(k - ii - 1) -= nu(ii);
        nu(ii) = 0.0;
      }
      if (nu(k - ii - 1) < 0.0) {
        nu(ii) -= nu(k - ii - 1);
        nu(k - ii - 1) = 0.0;
      }
    }
  }
  arma::mat levy_par = join_rows(basis, nu);
  
  return levy_par;
}
