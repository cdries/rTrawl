#include "RcppArmadillo.h"
#include "trawl_wrap.h"
#include "observe_process.h"
#include "cum.h"
#include "acf_ccf_helper.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export()]]
arma::vec acf_sample_p(double h, arma::vec x_grid, arma::vec p_grid, double TT, int lag_max) {
  // computes sample autocorrelation function of the process
  // 
  // arguments:
  // h        : time interval between observations / observation frequency
  // x_grid   : vector with times for given process values 
  // p_grid   : vector with process values
  // TT       : endpoint of observation interval
  // lag_max  : number of lags to look back, excludes lag 0
  //
  // author: Dries Cornilly
  
  arma::vec acfh = arma::zeros(lag_max);
  
  // centered observations
  p_grid -= cum_sample(1, x_grid, p_grid, TT);
  
  // variance
  double k2L = cum_sample(2, x_grid, p_grid, TT);
  
  // autocovariance
  for (int ii = 0; ii < lag_max; ii++) acfh(ii) = ccf_helper(x_grid, p_grid, x_grid, p_grid,TT, (ii + 1) * h);
  
  // as correlations
  acfh /= k2L;
  
  return acfh;
}

// [[Rcpp::export()]]
arma::vec acf_sample_dp(double h, arma::vec x_grid, arma::vec p_grid, 
                        double T0, double TT, int lag_max, int multi) {
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
  
  arma::vec T0_offset = arma::linspace(0.0, h, multi + 1);
  arma::vec acfh = arma::zeros(lag_max);
  double k2L;
  double k1L;
  for (int mm = 0; mm < multi; mm++) {
    
    // observe process
    List obs_process = observe_process(x_grid, p_grid, T0 + T0_offset(mm), TT, h);
    arma::vec obs_x = obs_process["x_grid_observed"];
    arma::vec obs_p = obs_process["p_grid_observed"];
    
    // take differences
    arma::vec diff_p = arma::diff(obs_p);

    // mean center the difference observations
    k1L = cum_dp_sample(1, h, obs_x, diff_p, TT);

    // variance
    k2L = cum_dp_sample(2, h, obs_x, diff_p, TT);

    // autocorrelations
    for (int ii = 0; ii < lag_max; ii++) {
      acfh(ii) += (ccf_dp_helper(h, obs_x, diff_p, obs_x, diff_p, TT, ii + 1) - k1L * k1L) / k2L;
    }
  }
  
  acfh /= multi;

  return acfh;
}

// [[Rcpp::export()]]
arma::vec acf_trawl_p(double h, std::string trawl, arma::vec trawl_par, int lag_max) {
  // computes theoretical autocorrelation function of the process
  // 
  // arguments:
  // h          : time interval between observations / observation frequency
  // trawl      : the trawl function in the process
  // trawl_par  : trawl parameters
  // lag_max    : number of lags to look back, excludes lag 0
  //
  // author: Dries Cornilly
  
  arma::vec h_vec = arma::linspace(0.0, lag_max * h, lag_max + 1);
  arma::vec acfh = leb_AtA(h_vec, trawl, trawl_par);
  acfh /= acfh(0);
  acfh = acfh.tail(lag_max);
  
  return acfh;
}

// [[Rcpp::export()]]
arma::vec acf_trawl_dp(double h, std::string trawl, 
                       arma::vec trawl_par,  double b, int lag_max) {
  // computes theoretical autocorrelation function of the process
  // 
  // arguments:
  // h          : length of differences / observation frequency
  // trawl      : the trawl function in the process
  // trawl_par  : trawl parameters
  // b          : b parameter as used in Shepard and Yang (2017)
  // lag_max    : number of lags to look back, excludes lag 0
  //
  // author: Dries Cornilly
  
  arma::vec h_vec = arma::linspace(0.0, (lag_max + 1.0) * h, lag_max + 2);
  arma::vec leb = leb_AtA(h_vec, trawl, trawl_par);
  
  arma::vec acfh(lag_max);
  if (b > 1.0 - std::numeric_limits<double>::epsilon()) {
    acfh = arma::zeros(lag_max);
  } else {
    acfh = (-leb.tail(lag_max) + 2.0 * leb.subvec(1, lag_max) - leb.head(lag_max)) / 
      (2.0 * (leb(0) - leb(1)) + h * b / (1.0 - b));
  }
  
  return acfh;
}

// [[Rcpp::export()]]
List acf_BN_V(double h, std::string trawl, arma::vec trawl_par, int lag_max) {
  // computes theoretical autocorrelation function of the process
  // together with the gradient with respect to the trawl parameters
  // 
  // arguments:
  // h          : time interval between observations / observation frequency
  // trawl      : the trawl function in the process
  // trawl_par  : trawl parameters
  // lag_max    : number of lags to look back, excludes lag 0
  //
  // author: Dries Cornilly
  
  // get theoretical acf
  arma::vec acf_theor = acf_trawl_p(h, trawl, trawl_par, lag_max);
  
  // compute necessary lebesgue measures and derivatives
  arma::vec h_vec = arma::linspace(0.0, lag_max * h, lag_max + 1);
  
  arma::mat dlebh = d_leb_AtA(h_vec, trawl, trawl_par);
  arma::mat dleb0 = arma::repmat(dlebh.row(0), lag_max, 1);
  
  arma::mat lebh = arma::repmat(leb_AtA(h_vec, trawl, trawl_par), 1, dlebh.n_cols);
  arma::mat leb0 = arma::repmat(lebh.row(0), lag_max, 1);
  arma::mat leb0_2inv = arma::pow(leb0, -2.0);
  
  // combine into the gradient
  arma::mat acf_grad = (dlebh.tail_rows(lag_max) % leb0 - lebh.tail_rows(lag_max) % dleb0) % leb0_2inv;
  
  List out;
  out["acf_theor"] = acf_theor;
  out["acf_grad"] = acf_grad;
  
  return out;
}
