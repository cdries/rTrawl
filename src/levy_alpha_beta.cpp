#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export()]]
List levy_alpha_beta(arma::vec p_grid, double T0, double TT) {
  // computes instantaneous jump probabilities and
  // jump rate, as in Shephard and Yang (2017) - needs integer jump sizes!
  // 
  // arguments:
  // p_grid   : vector with process values
  // T0       : beginpoint of observation interval
  // TT       : endpoint of observation interval
  //
  // author: Dries Cornilly
  
  double Mprec = std::numeric_limits<double>::epsilon();
  
  arma::vec Dp_grid = diff(p_grid, 1);
  double beta_0 = sum(abs(diff(p_grid)) > sqrt(Mprec)) / (TT - T0);
  
  // build basis based on unique observed jump sizes
  arma::vec LB = sort(unique(Dp_grid));
  int n_basis = LB.n_elem;
  arma::mat alpha_LB_temp(n_basis, 2);
  int tt = 0;
  for (int ii = 0; ii < n_basis; ii++) {
    if (LB(ii) * LB(ii) > 0.0001) {
      alpha_LB_temp(tt, 0) = LB(ii);
      alpha_LB_temp(tt, 1) = sum((Dp_grid > LB(ii) - 0.5) && (Dp_grid < LB(ii) + 0.5));
      tt++;
    }
  }
  arma::mat alpha_LB = alpha_LB_temp.submat(0, 0, tt - 1, 1);
  alpha_LB.col(1) /= sum(alpha_LB.col(1));
  
  List out;
  out["beta_0"] = beta_0;
  out["levy_alpha"] = alpha_LB;
  
  return out;
}
