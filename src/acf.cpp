#include "RcppArmadillo.h"
#include "trawl_wrap.h"
#include "observe_process.h"
#include "cum.h"
#include "acf_ccf_helper.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export()]]
arma::vec acf_sample_p(double h, arma::vec x_grid, arma::vec p_grid, double TT, int lag_max) {
  
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
  
  arma::vec T0_offset = arma::linspace(0.0, h, multi + 1);
  arma::vec acfh = arma::zeros(lag_max);
  
  for (int mm = 0; mm < multi; mm++) {
    
    // observe process
    List obs_process = observe_process(x_grid, p_grid, T0 + T0_offset(mm), TT, h);
    arma::vec obs_x = obs_process["x_grid_observed"];
    arma::vec obs_p = obs_process["p_grid_observed"];
    
    // observation grid of the returns
    int n = obs_x.n_elem;
    obs_x = obs_x.tail(n - 1);
    n--;
    
    // mean and centered observations
    arma::vec diff_p = arma::diff(obs_p);
    double n_full = floor((TT - T0 - T0_offset(mm)) / h);
    double k1 = arma::sum(diff_p) / n_full;
    diff_p -= k1;
    
    // variance
    double k2 = arma::sum(diff_p % diff_p) / n_full;
    
    // autocorrelation
    for (int ii = 1; ii < lag_max + 1; ii++) {
      double ac = 0.0;
      for (int tt = 1; tt < n; tt++) {
        int kk = tt - 1;
        while ((kk > -0.5) && (tt - kk < ii + 0.5)) {
          if ((obs_x(kk) > obs_x(tt) - (ii + 1.0) * h) && (obs_x(kk) < obs_x(tt) - (ii - 1.0) * h)) {
            ac += obs_p(kk) * obs_x(tt);
            kk = -1;
          } else {
            kk++;
          }
        }
      }
      ac /= n_full - ii;
      acfh(ii) += ac / k2;
    }
  }
  acfh /= multi;
  
  return acfh;
}

// [[Rcpp::export()]]
arma::vec acf_trawl_p(double h, std::string trawl, arma::vec trawl_par, int lag_max) {
  
  arma::vec h_vec = arma::linspace(0.0, lag_max * h, lag_max + 1);
  arma::vec acfh = leb_AtA(h_vec, trawl, trawl_par);
  acfh /= acfh(0);
  acfh = acfh.tail(lag_max);
  
  return acfh;
}

// [[Rcpp::export()]]
arma::vec acf_trawl_dp(double h, std::string trawl, 
                       arma::vec trawl_par,  double b, int lag_max) {
  
  arma::vec h_vec = arma::linspace(0.0, (lag_max + 1.0) * h, lag_max + 2);
  arma::vec leb = leb_AtA(h_vec, trawl, trawl_par);
  
  arma::vec acfh(lag_max);
  if (b > 1.0 - std::numeric_limits<double>::epsilon()) {
    acfh = arma::zeros(lag_max);
  } else {
    acfh = (-leb.tail(lag_max) + 2.0 * leb.subvec(1, lag_max) - leb.head(lag_max)) / 
      (2.0 * (leb(1) - leb(0)) + h * b / (1.0 - b));
  }
  
  return acfh;
}

// [[Rcpp::export()]]
List acf_BN_V(double h, std::string trawl, arma::vec trawl_par, int lag_max) {
  
  arma::vec acf_theor = acf_trawl_p(h, trawl, trawl_par, lag_max);
  
  arma::vec h_vec = arma::linspace(0.0, lag_max * h, lag_max + 1);
  
  arma::mat dlebh = d_leb_AtA(h_vec, trawl, trawl_par);
  arma::mat dleb0 = arma::repmat(dlebh.row(0), lag_max, 1);
  
  arma::mat lebh = arma::repmat(leb_AtA(h_vec, trawl, trawl_par), 1, dlebh.n_cols);
  arma::mat leb0 = arma::repmat(lebh.row(0), lag_max, 1);
  arma::mat leb0_2inv = arma::pow(leb0, -2.0);
  
  arma::mat acf_grad = (dlebh.tail_rows(lag_max) % leb0 - lebh.tail_rows(lag_max) % dleb0) % leb0_2inv;
  
  List out;
  out["acf_theor"] = acf_theor;
  out["acf_grad"] = acf_grad;
  
  return out;
}
