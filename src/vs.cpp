#include "RcppArmadillo.h"
#include "trawl_wrap.h"
#include "observe_process.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export()]]
arma::vec vs_sample(arma::vec h, arma::vec x_grid, arma::vec p_grid, double T0, double TT, int multi) {
  // sample variance signature values, see Shephard and Yang (2017)
  //
  // arguments:
  // h        : length of differences / observation frequency
  // x_grid   : vector with times for given process values 
  // p_grid   : vector with process values
  // T0       : beginpoint of observation interval
  // TT       : endpoint of observation interval
  // multi    : number of offsets to use in order to reduce estimation variance
  //
  // author: Dries Cornilly
  
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

// [[Rcpp::export()]]
List vs_C(arma::vec h, std::string trawl, arma::vec trawl_par, double omega,
          double xi, bool include_b, double eta, bool include_cum1) {
  // theoretical variance signature plot and its gradient, see XXX
  //
  // arguments:
  // h          : length of differences / observation frequency
  // trawl      : trawl of the process
  // trawl_par  : values for the trawl
  // omega      : second cumulant of the noise process
  // xi         : second cumulant of the efficient process
  // include_b  : whether or not to include an efficient part (xi)
  // eta        : squared first cumulant of the efficient part
  // include_cum1 : whether or not to include the first cumulant 
  //
  // author: Dries Cornilly
  
  int n_h = h.n_elem;
  
  arma::vec hinv = arma::pow(h, -1.0);
  arma::vec leb0 = arma::ones(n_h) * leb_AtA(arma::zeros(1), trawl, trawl_par)(0);
  arma::vec lebh = leb_AtA(h, trawl, trawl_par);
  arma::vec vs_theor = 2.0 * (leb0 - lebh) % hinv * omega + xi + h * eta;
  
  arma::mat hinv_mat = arma::repmat(hinv, 1, number_parameters_trawl(trawl));
  arma::mat dleb0 = arma::repmat(d_leb_AtA(arma::zeros(1), trawl, trawl_par), n_h, 1);
  arma::mat dlebh = d_leb_AtA(h, trawl, trawl_par);
  arma::mat vs_grad = 2.0 * (dleb0 - dlebh) % hinv_mat * omega;
  vs_grad = arma::join_rows(vs_grad, 2.0 * (leb0 - lebh) % hinv);
  if (include_b) vs_grad = arma::join_rows(vs_grad, arma::ones(n_h));
  if (include_cum1) vs_grad = arma::join_rows(vs_grad, h);
  
  List out;
  out["vs_theor"] = vs_theor;
  out["vs_grad"] = vs_grad;
  
  return out;
}

// [[Rcpp::export()]]
List vs_SY(arma::vec h, std::string trawl, arma::vec trawl_par, double beta_0,
           arma::mat levy_alpha, bool include_cum1, double b, bool include_b) {
  // theoretical variance signature plot and its gradient, see Shephard and Yang (2017)
  //
  // arguments:
  // h          : length of differences / observation frequency
  // trawl      : trawl of the process
  // trawl_par  : values for the trawl
  // beta_0     : rate of observed price changes, see Shephard and Yang (2017)
  // levy_alpha : instantaneous jump probabilities, see Shephard and Yang (2017)
  // include_cum1 : whether or not to include the first cumulant
  // b          : b value used, see Shephard and Yang (2017)
  // include_b  : whether or not to include an efficient part (xi)
  //
  // author: Dries Cornilly
  
  int n_h = h.n_elem;
  
  double k1Lpart = 0.0;
  if (include_cum1) k1Lpart = beta_0 * arma::sum(levy_alpha.col(0) % levy_alpha.col(1));
  double k2L = beta_0 * arma::sum(levy_alpha.col(0) % levy_alpha.col(0) %
                                  levy_alpha.col(1)) / (2.0 - b);
  arma::vec hinv = arma::pow(h, -1.0);
  arma::vec leb0 = arma::ones(n_h) * leb_AtA(arma::zeros(1), trawl, trawl_par)(0);
  arma::vec lebh = leb_AtA(h, trawl, trawl_par);
  arma::vec vs_theor = (b + 2.0 * (1.0 - b) * (leb0 - lebh) % hinv) * k2L + h * k1Lpart * k1Lpart;
  
  arma::mat hinv_mat = arma::repmat(hinv, 1, number_parameters_trawl(trawl));
  arma::mat dleb0 = arma::repmat(d_leb_AtA(arma::zeros(1), trawl, trawl_par), n_h, 1);
  arma::mat dlebh = d_leb_AtA(h, trawl, trawl_par);
  arma::mat vs_grad = 2.0 * (1.0 - b) * (dleb0 - dlebh) % hinv_mat * k2L;
  if (include_b) {
    vs_grad = arma::join_rows(vs_grad,  vs_theor / (2.0 - b) +
      (1.0 - 2.0 * (leb0 - lebh) % hinv) * k2L);
  }
  
  List out;
  out["vs_theor"] = vs_theor;
  out["vs_grad"] = vs_grad;
  
  return out;
}
