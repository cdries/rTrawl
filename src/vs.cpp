#include "RcppArmadillo.h"
#include "trawl_wrap.h"
#include "observe_process.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export()]]
arma::vec vs_sample(arma::vec h, arma::vec x_grid, arma::vec p_grid, double T0, double TT, int multi) {
  
  int n_h = h.n_elem;
  double DTT = TT - T0;
  arma::vec vs = arma::zeros(n_h);
  
  for (int ii = 0; ii < n_h; ii++) {
    arma::vec T0_offset = arma::linspace(0.0, 0.5, multi + 1);
    
    double vs_ii = 0.0;
    for (int mm = 0; mm < multi; mm++) {
      arma::vec p_grid_h = observe_process(x_grid, p_grid, T0 + T0_offset(mm), TT, h(ii))["p_grid_observed"];
      arma::vec Dp_grid_h = arma::diff(p_grid_h);
      vs_ii += arma::sum(Dp_grid_h % Dp_grid_h);
    }
    vs(ii) = vs_ii / (multi * DTT);
  }
  
  return vs;
}

// // [[Rcpp::export()]]
// List vs_SY(arma::vec h, std::string trawl, arma::vec trawl_par, double beta_0, 
//            arma::mat levy_alpha, bool include_cum1, double b, bool include_b) {
//   
//   int n_h = h.n_elem;
//   
//   double k1Lpart = 0.0;
//   if (include_cum1) k1Lpart = beta_0 * arma::sum(levy_alpha.col(0) % levy_alpha.col(1));
//   double k2L = beta_0 * arma::sum(levy_alpha.col(0) % levy_alpha.col(0) % 
//                                   levy_alpha.col(1)) / (2.0 - b);
//   arma::vec hinv = arma::pow(h, -1.0);
//   arma::vec leb0 = arma::ones(n_h) * leb_AtA(0.0, trawl, trawl_par, b)(0);
//   arma::vec vs_theor = (b + 2.0 * (leb0 - leb_AtA(h, trawl, trawl_par, b)) * hinv) * k2L + h * k1Lpart * k1Lpart;
//   
//   List out;
//   out["vs_theor"] = vs_theor;
//   out["vs_grad"] = 0.0;
//   
//   return out;
// }
