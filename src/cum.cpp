#include "RcppArmadillo.h"
#include "trawl_wrap.h"
#include "acf_ccf_helper.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export()]]
double cum_sample(int ord, arma::vec x_grid, arma::vec p_grid, double TT) {
  // computes sample cumulants of an observed process
  // 
  // arguments:
  // ord      : order of the moment - 1 or 2
  // x_grid   : vector with times for given process values
  // p_grid   : vector with process values
  // TT       : endpoint of observation interval
  //
  // author: Dries Cornilly
  
  // mean
  int n = x_grid.n_elem;
  arma::vec diff_x = arma::diff(x_grid);
  double k1L = (arma::sum(diff_x % p_grid.head(n - 1)) +
                (TT - x_grid(n - 1)) * p_grid(n - 1)) / (TT - x_grid(0));
  
  double cum = 0.0;
  if (ord == 1) {
    cum = k1L;
    
  } else if (ord == 2) {
    p_grid -= k1L;
    
    // variance
    cum = (arma::sum(diff_x % p_grid.head(n - 1) % p_grid.head(n - 1)) +
      (TT - x_grid(n - 1)) * p_grid(n - 1) * p_grid(n - 1)) / (TT - x_grid(0));
  }
  
  return cum;
}

double cum_dp_sample(int ord, double h, arma::vec x_grid, arma::vec diff_p, double TT) {
  // computes sample cumulants of a differenced process
  // 
  // arguments:
  // ord      : order of the moment - 1 or 2
  // h        : length of differences / observation frequency
  // x_grid   : vector with times for given process values
  // diff_p   : vector with observed process differences
  // TT       : endpoint of observation interval
  //
  // author: Dries Cornilly
  
  int n = diff_p.n_elem;
  double n_total = floor((TT - x_grid(0)) / h);
  x_grid = x_grid.tail(n);
  double k1L = arma::sum(diff_p) / n_total;
  
  double cum = 0.0;
  if (ord == 1) {
    cum = k1L;
  } else if (ord == 2) {
    cum = arma::sum(diff_p % diff_p) / n_total - k1L * k1L;
  }
  
  return cum;
}

// [[Rcpp::export()]]
List levy_cum_mv2fit(double T0, double TT, List x_grid, List p_grid, List trawl, List trawl_par, int p) {
  // computes univariate first cumulants and multivariate second cumulants of 
  // the underlying multivariate levy seed - NOT the sample moments
  // 
  // arguments:
  // T0         : beginpoint of observation interval
  // TT         : endpoint of observation interval
  // x_grid     : List where each vector contains times for given process values
  // p_grid     : List where each vector contains process values
  // trawl      : List with trawl for each component
  // trawl_par  : List with trawl parameters for each component
  // p          : dimension of the process
  //
  // author: Dries Cornilly
  
  // initialize
  arma::vec k1_sample = arma::zeros(p);
  arma::vec k2_sample = arma::zeros(p * (p + 1) / 2);
  
  int kk = 0;
  for (int ii = 0; ii < p; ii++) {
    
    // extract process1 from list
    std::string trawl1 = trawl(ii);
    arma::vec trawl_par1 = trawl_par(ii);
    arma::vec p_grid1 = p_grid(ii);
    arma::vec x_grid1 = x_grid(ii);
    
    // compute mean and standardize
    arma::vec lebAii = leb_AtA(arma::zeros(1), trawl1, trawl_par1);
    double k1_temp = cum_sample(1, x_grid1, p_grid1, TT);
    k1_sample(ii) = k1_temp / lebAii(0);
    p_grid1 -= k1_temp;
    
    for (int jj = ii; jj < p; jj++) {
      
      // extract process 2 from list
      arma::vec p_grid2 = p_grid(jj);
      arma::vec x_grid2 = x_grid(jj);
      p_grid2 -= cum_sample(1, x_grid2, p_grid2, TT);
      
      // compute second cumulant, both univariate and multivariate
      k2_sample(kk) = ccf_helper(x_grid1, p_grid1, x_grid2, p_grid2, TT, 0.0) /
        leb_autocorrelator(arma::zeros(1), trawl1, trawl_par1, trawl(jj), trawl_par(jj))(0);
      kk++;
    }
  }
  
  List out;
  out["k1_sample"] = k1_sample;
  out["k2_sample"] = k2_sample;

  return out;
}
