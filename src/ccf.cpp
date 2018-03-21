#include "RcppArmadillo.h"
#include "trawl_wrap.h"
#include "observe_process.h"
#include "cum.h"
#include "acf_ccf_helper.h"
#include "levy_wrap.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export()]]
arma::vec ccf_sample_p(double h, arma::vec x_grid1, arma::vec p_grid1, arma::vec x_grid2, 
                       arma::vec p_grid2, double TT, int lag_max) {
  
  int n_lags = 2 * lag_max + 1;
  arma::vec ccfh = arma::zeros(n_lags);
  
  // centered observations
  p_grid1 -= cum_sample(1, x_grid1, p_grid1, TT);
  p_grid2 -= cum_sample(1, x_grid2, p_grid2, TT);
  
  // standard deviations
  double sd1 = sqrt(cum_sample(2, x_grid1, p_grid1, TT));
  double sd2 = sqrt(cum_sample(2, x_grid2, p_grid2, TT));
  
  // cross covariance
  for (int ii = 0; ii < lag_max; ii++) {
    ccfh(ii) = ccf_helper(x_grid1, p_grid1, x_grid2, p_grid2, TT, (lag_max - ii) * h);
  }
  for (int ii = lag_max; ii < n_lags; ii++) {
    ccfh(ii) = ccf_helper(x_grid2, p_grid2, x_grid1, p_grid1, TT, (ii - lag_max) * h);
  }
  
  // as correlations
  ccfh /= sd1 * sd2;
  
  return ccfh;
}

// [[Rcpp::export()]]
arma::vec ccf_sample_dp(double h, arma::vec x_grid1, arma::vec p_grid1, arma::vec x_grid2, 
                        arma::vec p_grid2, double T0, double TT, int lag_max, int multi) {
  
  arma::vec T0_offset = arma::linspace(0.0, h, multi + 1);
  int n_lags = 2 * lag_max + 1;
  arma::vec ccfh = arma::zeros(n_lags);
  
  for (int mm = 0; mm < multi; mm++) {
    
    // observe process
    List obs_process = observe_process(x_grid1, p_grid1, T0 + T0_offset(mm), TT, h);
    arma::vec obs_x1 = obs_process["x_grid_observed"];
    arma::vec obs_p1 = obs_process["p_grid_observed"];
    obs_process = observe_process(x_grid2, p_grid2, T0 + T0_offset(mm), TT, h);
    arma::vec obs_x2 = obs_process["x_grid_observed"];
    arma::vec obs_p2 = obs_process["p_grid_observed"];
    
    // take differences
    arma::vec diff_p1 = arma::diff(obs_p1);
    arma::vec diff_p2 = arma::diff(obs_p2);
    
    // mean center the differences observations
    double mu1 = cum_dp_sample(1, h, obs_x1, diff_p1, TT);
    double mu2 = cum_dp_sample(1, h, obs_x2, diff_p2, TT);
    
    // standard deviations
    double sd1 = sqrt(cum_dp_sample(2, h, obs_x1, diff_p1, TT));
    double sd2 = sqrt(cum_dp_sample(2, h, obs_x2, diff_p2, TT));
    
    // cross correlations
    for (int ii = 0; ii < lag_max; ii++) {
      ccfh(ii) += (ccf_dp_helper(h, obs_x1, diff_p1, obs_x2, diff_p2, TT, lag_max - ii) - mu1 * mu2) / (sd1 * sd2);
    }
    for (int ii = lag_max; ii < n_lags; ii++) {
      ccfh(ii) += (ccf_dp_helper(h, obs_x2, diff_p2, obs_x1, diff_p1, TT, ii - lag_max) - mu1 * mu2) / (sd1 * sd2);
    }
  }
  
  ccfh /= multi;
  
  return ccfh;
}

// [[Rcpp::export()]]
arma::vec ccf_trawl_p(double h, std::string trawl1, arma::vec trawl_par1, 
                      std::string trawl2, arma::vec trawl_par2, 
                      std::string levy_seed, arma::mat levy_par, arma::mat design_matrix,
                      int lag_max) {
  
  arma::vec h_vec = arma::linspace(-h * lag_max, h * lag_max, 2 * lag_max + 1);
  arma::vec ccfh = leb_autocorrelator(h_vec, trawl1, trawl_par1, trawl2, trawl_par2);
  
  arma::vec leb1 = leb_AtA(arma::zeros(1), trawl1, trawl_par1);
  arma::vec leb2 = leb_AtA(arma::zeros(1), trawl2, trawl_par2);
  
  arma::mat varcovar = levy_varcovar(levy_seed, levy_par, design_matrix);
  
  ccfh *= varcovar(0, 1) / sqrt(leb1 * leb2 * varcovar(1, 1) * varcovar(0, 0));
  
  return ccfh;
}

// [[Rcpp::export()]]
arma::vec ccf_trawl_dp(double h, std::string trawl1, arma::vec trawl_par1, 
                       std::string trawl2, arma::vec trawl_par2, arma::vec b,
                       std::string levy_seed, arma::mat levy_par, arma::mat design_matrix, 
                       int lag_max) {
  
  arma::vec ccfh = arma::zeros(5);
  
  return ccfh;
}
