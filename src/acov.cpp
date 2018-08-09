#include "RcppArmadillo.h"
#include "observe_process.h"
#include "cum.h"
#include "acf_ccf_helper.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export()]]
List acov(double h, arma::vec x_grid, arma::vec p_grid, 
               double T0, double TT, int lag_max) {
  // computes sample autocorrelation function of the differenced process
  // 
  // arguments:
  // h        : length of differences / observation frequency
  // x_grid   : vector with times for given process values 
  // p_grid   : vector with process values
  // T0       : beginpoint of observation interval
  // TT       : endpoint of observation interval
  // lag_max  : number of lags to look back, excludes lag 0
  // multi    : number of offsets to use in order to reduce estimation variance
  //
  // author: Dries Cornilly
  
  arma::vec acv = arma::zeros(lag_max + 1);
  double k2L;
  double k1L;
  double DTT = (TT - T0) * (TT - T0);
  
  // observe process
  List obs_process = observe_process(x_grid, p_grid, T0, TT, h);
  arma::vec obs_x = obs_process["x_grid_observed"];
  arma::vec obs_p = obs_process["p_grid_observed"];
  
  double n_total = floor((TT - x_grid(0)) / h);
  
  // take differences
  arma::vec diff_p = arma::diff(obs_p);
  arma::vec diff_p2 = diff_p % diff_p;
  
  // mean center the difference observations
  k1L = cum_dp_sample(1, h, obs_x, diff_p2, TT);
  
  // variance
  k2L = cum_dp_sample(2, h, obs_x, diff_p2, TT);
  
  // build covariance
  acv(0) += n_total * k2L;
  for (int ii = 0; ii < lag_max; ii++) {
    acv(ii + 1) += 2.0 * (n_total - ii - 1) * 
      (ccf_dp_helper(h, obs_x, diff_p2, obs_x, diff_p2, TT, ii + 1) - k1L * k1L);
  }
  
  acv /= DTT;
  
  List out;
  out["acv"] = acv;

  return out;
}