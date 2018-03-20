#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// computes cross covariance between irregularly observed time series / h positive
// this means time series 2 is trailing: int p_1(t) p_2(t - h) dt
double ccf_helper(arma::vec x1_grid, arma::vec p1_grid, arma::vec x2_grid, 
                  arma::vec p2_grid, double TT, double h) {
  
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

// computes cross product between observed time series on h grid / h positive
// with lg a positive (or zero) integer
// this means time series 2 is trailing: sum p_1(t) p_2(t - lg * h)
double ccf_dp_helper(double h, arma::vec x_grid1, arma::vec diff_p1, arma::vec x_grid2, 
                     arma::vec diff_p2, double TT, int lg) {
  
  int n1 = diff_p1.n_elem;
  int n2 = diff_p2.n_elem;
  double n_total = floor((TT - x_grid2(0)) / h) - lg;
  arma::vec x1_grid = x_grid1.tail(n1);
  arma::vec x2_grid = x_grid2.tail(n2);
  
  // initialize
  double cc = 0.0;
  int ii1 = 0;
  int ii2 = 0;
  
  // iterate over time indices
  while (ii1 < n1 && ii2 < n2) {
     if (x2_grid(ii2) + h * lg < x1_grid(ii1) - 0.5 * h) {
       ii2++;
     } else  {
       if ((x2_grid(ii2) + h * lg < x1_grid(ii1) + 0.5 * h) && (x2_grid(ii2) + h * lg > x1_grid(ii1) - 0.5 * h)) {
         cc += diff_p1(ii1) * diff_p2(ii2);
       }
       ii1++;
     }
  }
  
  // divide over number of observations
  cc /= n_total;
  
  return cc;
}

