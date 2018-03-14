#include "RcppArmadillo.h"
#include "trawl_wrap.h"
#include "observe_process.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


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
