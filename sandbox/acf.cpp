#include "RcppArmadillo.h"
#include "trawl_wrap.h"
#include "observe_process.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


double acf_helper(arma::vec x_grid, arma::vec p_grid, double h, int lag) {
  
  double ac = 0.0;
  int n = x_grid.n_elem;
  
  // iterate over observations
  for (int ii = 1; ii < n; ii++) {
    
    // update at the new observation
    int Tcurrent = x_grid(ii);
    int kk = 1;
    while ((ii - kk > 0.5) && (x_grid(ii - kk) + lag * h > Tcurrent + 0.5 * h)) kk++;
    if (x_grid(ii - kk) + lag * h < Tcurrent + 0.5 * h) {
      
      // add most recent autocorrelation - where price has just changed
      ac += p_grid(ii - kk) * p_grid(ii); 
      
      // update autocorrelations strictly between t_ii and t_{ii = 1}
      while (Tcurrent > x_grid(ii - 1) + 1.5 * h) {
        if (x_grid(ii - kk) + lag * h < x_grid(ii - 1) + 1.5 * h) {
          ac += (Tcurrent - x_grid(ii - 1) - h) * p_grid(ii - kk) * p_grid(ii - 1);
          Tcurrent = x_grid(ii - 1);
        } else {
          ac += (Tcurrent - (x_grid(ii - kk) + h * lag)) * p_grid(ii - kk) * p_grid(ii - 1);
          Tcurrent = x_grid(ii - kk) + lag * h;
          kk++;
          if (ii - kk < -0.5) Tcurrent = x_grid(ii - 1);
        }
      }
    }
  }
  ac /= (x_grid(n - 1) - x_grid(0) - h * lag);
  
  return ac;
}


// [[Rcpp::export()]]
arma::vec acf_sample_p(double h, arma::vec x_grid, arma::vec p_grid, 
                       double T0, double TT, int lag_max, int multi) {
  
  arma::vec T0_offset = arma::linspace(0.0, h, multi + 1);
  arma::vec acfh = arma::zeros(lag_max);
  
  for (int mm = 0; mm < multi; mm++) {

    // observe process
    List obs_process = observe_process(x_grid, p_grid, T0 + T0_offset(mm), TT, h);
    arma::vec obs_x = obs_process["x_grid_observed"];
    arma::vec obs_p = obs_process["p_grid_observed"];

    // mean and centered observations
    int n = obs_x.n_elem;
    arma::vec diff_x = arma::diff(obs_x);
    double k1L = (arma::sum(diff_x % obs_p.head(n - 1)) +
                  (TT - obs_x(n - 1)) * obs_p(n - 1)) / (TT - obs_p(0));
    obs_p -= k1L;

    // variance
    double k2L = (arma::sum(diff_x % obs_p.head(n - 1) % obs_p.head(n - 1)) +
                  (TT - obs_x(n - 1)) * obs_p(n - 1) * obs_p(n - 1)) / (TT - obs_p(0));

    // autocorrelation
    for (int ii = 0; ii < lag_max; ii++) acfh(ii) += acf_helper(obs_x, obs_p, h, ii + 1) / k2L;
  }
  acfh /= multi;
  
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
