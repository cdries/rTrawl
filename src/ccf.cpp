#include "RcppArmadillo.h"
#include "trawl_wrap.h"
#include "observe_process.h"
#include "cum.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// computes cross covariance between irregularly observed time series / h positive
// this means time series 2 is trailing: int p_1(t) p_2(t - h) dt
double ccf_helper(arma::vec x1_grid, arma::vec p1_grid, arma::vec x2_grid, arma::vec p2_grid, 
                  double TT, double h) {
  
  int n1 = x1_grid.n_elem;
  int n2 = x2_grid.n_elem;
  
  // attach TT to time series 1 and 2 if required
  if (x1_grid(n1 - 1) < TT - std::numeric_limits<double>::epsilon()) {
    x1_grid = arma::join_cols(x1_grid, arma::ones(1) * TT);
    p1_grid = arma::join_cols(p1_grid, arma::ones(1) * p1_grid(n1 - 1));
    n1++;
  }
  if (x2_grid(n2 - 1) < TT - std::numeric_limits<double>::epsilon()) {
    x2_grid = arma::join_cols(x2_grid, arma::ones(1) * TT);
    p2_grid = arma::join_cols(p2_grid, arma::ones(1) * p2_grid(n2 - 1));
    n2++;
  }
  
  // initialize
  double cc = 0.0;
  double totalT = 0.0;
  int ii1 = 0;
  int ii2 = 0;
  double x1_current = x1_grid(ii1);
  double x2_current = x2_grid(ii1);
  double p1_current = p1_grid(ii2);
  double p2_current = p2_grid(ii2);
  
  // remove unnecessary time series 2 observations
  while (x2_grid(ii2) + h < x1_grid(ii1)) ii2++;
  if (ii2 > 0) {
    ii2--;
    x2_current = x1_grid(ii1) - h;
    p2_current = p2_grid(ii2);
  } else {
    
    // move time series 1 forward to 1st useful time and price
    while (x1_grid(ii1) < x2_grid(0) + h) ii1++;
    if (ii1 > 0) {
      ii1--;
      x1_current = x2_grid(ii2) + h;
      x2_current = x2_grid(ii2);
      p1_current = p1_grid(ii1);
      p2_current = p2_grid(ii2);
    }
  }
  
  // iterate over time indices
  while (ii1 < n1 - 1 && ii2 < n2 - 1) {

    if (x1_grid(ii1 + 1) > x2_grid(ii2 + 1) + h) {

      // jump in time1 is too large to make at once
      double deltaT = (x2_grid(ii2 + 1) - x2_current);
      cc += deltaT * p1_current * p2_current;
      ii2++;
      x2_current = x2_grid(ii2);
      p2_current = p2_grid(ii2);
      x1_current += deltaT;
      totalT += deltaT;
    } else {

      // jump in time2 is too large to make at once
      double deltaT = (x1_grid(ii1 + 1) - x1_current);
      cc += deltaT * p1_current * p2_current;
      ii1++;
      x1_current = x1_grid(ii1);
      p1_current = p1_grid(ii1);
      x2_current += deltaT;
      totalT += deltaT;
    }
  }

  // divide over information time period
  cc /= totalT;

  return cc;
}


// [[Rcpp::export()]]
arma::vec ccf_sample_p(double h, arma::vec x_grid1, arma::vec p_grid1, arma::vec x_grid2, 
                       arma::vec p_grid2, double TT, int lag_max) {
  
  int n_lags = 2 * lag_max + 1;
  arma::vec ccfh = arma::zeros(n_lags);
  arma::vec h_vec = arma::linspace(-h * lag_max, h * lag_max, n_lags);
  
  // mean and centered observations
  p_grid1 -= cum_sample(1, x_grid1, p_grid1, TT);
  p_grid1 -= cum_sample(1, x_grid2, p_grid2, TT);
  
  // standard deviations
  double sd1 = sqrt(cum_sample(2, x_grid1, p_grid1, TT));
  double sd2 = sqrt(cum_sample(2, x_grid2, p_grid2, TT));
  
  // cross covariance
  for (int ii = 0; ii < lag_max; ii++) {
    ccfh(ii) = ccf_helper(x_grid1, p_grid1, x_grid2, p_grid2, TT, (lag_max - ii) * h);
  }
  for (int ii = lag_max; ii < n_lags; ii++) {
    ccfh(ii) = ccf_helper(x_grid2, p_grid2, x_grid1, p_grid1, TT, (lag_max - ii) * h);
  }
  
  // as correlations
  ccfh /= sd1 * sd2;
  
  return ccfh;
}

// // [[Rcpp::export()]]
// arma::vec acf_sample_dp(double h, arma::vec x_grid, arma::vec p_grid, 
//                         double T0, double TT, int lag_max, int multi) {
//   
//   arma::vec T0_offset = arma::linspace(0.0, h, multi + 1);
//   arma::vec acfh = arma::zeros(lag_max);
//   
//   for (int mm = 0; mm < multi; mm++) {
//     
//     // observe process
//     List obs_process = observe_process(x_grid, p_grid, T0 + T0_offset(mm), TT, h);
//     arma::vec obs_x = obs_process["x_grid_observed"];
//     arma::vec obs_p = obs_process["p_grid_observed"];
//     
//     // observation grid of the returns
//     int n = obs_x.n_elem;
//     obs_x = obs_x.tail(n - 1);
//     n--;
//     
//     // mean and centered observations
//     arma::vec diff_p = arma::diff(obs_p);
//     double n_full = floor((TT - T0 - T0_offset(mm)) / h);
//     double k1 = arma::sum(diff_p) / n_full;
//     diff_p -= k1;
//     
//     // variance
//     double k2 = arma::sum(diff_p % diff_p) / n_full;
//     
//     // autocorrelation
//     for (int ii = 1; ii < lag_max + 1; ii++) {
//       double ac = 0.0;
//       for (int tt = 1; tt < n; tt++) {
//         int kk = tt - 1;
//         while ((kk > -0.5) && (tt - kk < ii + 0.5)) {
//           if ((obs_x(kk) > obs_x(tt) - (ii + 1.0) * h) && (obs_x(kk) < obs_x(tt) - (ii - 1.0) * h)) {
//             ac += obs_p(kk) * obs_x(tt);
//             kk = -1;
//           } else {
//             kk++;
//           }
//         }
//       }
//       ac /= n_full - ii;
//       acfh(ii) += ac / k2;
//     }
//   }
//   acfh /= multi;
// 
//   return acfh;
// }
// 
// // [[Rcpp::export()]]
// arma::vec acf_trawl_p(double h, std::string trawl, arma::vec trawl_par, int lag_max) {
//   
//   arma::vec h_vec = arma::linspace(0.0, lag_max * h, lag_max + 1);
//   arma::vec acfh = leb_AtA(h_vec, trawl, trawl_par);
//   acfh /= acfh(0);
//   acfh = acfh.tail(lag_max);
//   
//   return acfh;
// }
// 
// // [[Rcpp::export()]]
// arma::vec acf_trawl_dp(double h, std::string trawl, 
//                        arma::vec trawl_par,  double b, int lag_max) {
//   
//   arma::vec h_vec = arma::linspace(0.0, (lag_max + 1.0) * h, lag_max + 2);
//   arma::vec leb = leb_AtA(h_vec, trawl, trawl_par);
//   
//   arma::vec acfh(lag_max);
//   if (b > 1.0 - std::numeric_limits<double>::epsilon()) {
//     acfh = arma::zeros(lag_max);
//   } else {
//     acfh = (-leb.tail(lag_max) + 2.0 * leb.subvec(1, lag_max) - leb.head(lag_max)) / 
//       (2.0 * (leb(1) - leb(0)) + h * b / (1.0 - b));
//   }
//   
//   return acfh;
// }
