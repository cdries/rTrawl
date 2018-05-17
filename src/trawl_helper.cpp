#include "RcppArmadillo.h"
#include "trawl_wrap.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


double leb_EXP_EXP(double h, double lambda1, double lambda2) {
  // compute lebesgue measure of the minimum of trawl 1 and trawl 2 for exponential trawls,
  // h always positive! - trawl 1 is on (-infty, 0) - trawl 2 is on (-infty, h)
  //
  // arguments:
  // h          : argument for the support of trawl 2
  // lambda1    : trawl parameter for trawl 1
  // lambda2    : trawl parameter for trawl 2
  //
  // author: Dries Cornilly
  
  double leb;
  
  if (lambda2 < lambda1) {
    leb = exp(-lambda2 * h) / lambda2 - 
      (1 / lambda2 - 1 / lambda1) * exp(lambda1 * lambda2 * h / (lambda2 - lambda1));
  } else {
    leb = exp(-lambda2 * h) / lambda2;
  }
  
  return leb; 
}

// h always positive! - trawl 1 is on (-infty, 0) - trawl 2 is on (-infty, h)
double leb_GEN_GEN(double h, std::string trawl1, arma::vec trawl_par1, 
                   std::string trawl2, arma::vec trawl_par2) {
  // approximate lebesgue measure of the minimum of trawl 1 and trawl 2 for general trawls,
  // h always positive! - trawl 1 is on (-infty, 0) - trawl 2 is on (-infty, h)
  //
  // arguments:
  // h          : argument for the support of trawl 2
  // trawl1     : trawl 1
  // trawl_par1 : trawl parameters for trawl 1
  // trawl2     : trawl 2
  // trawl_par2 : trawl parameters for trawl 2
  //
  // author: Dries Cornilly
  
  double ub = 100.0;
  double Nsupp1 = 10000;
  double Nsupp2 = 100000;
  
  arma::vec supp0 = (arma::logspace(0.0, 1.0, Nsupp1) - 1.0) / (10.0 - 1.0);
  supp0 = arma::join_cols(arma::zeros(1), supp0);
  arma::vec supp1 = arma::linspace(1.0, ub, Nsupp2);
  
  arma::vec val0_trawl1 = trawl_function(-supp0, trawl1, trawl_par1);
  arma::vec val1_trawl1 = trawl_function(-supp1, trawl1, trawl_par1);
  
  arma::vec val0_trawl2 = trawl_function(-supp0 - h, trawl2, trawl_par2);
  arma::vec val1_trawl2 = trawl_function(-supp1 - h, trawl2, trawl_par2);
  
  arma::vec val0 = arma::min(val0_trawl1, val0_trawl2);
  arma::vec val1 = arma::min(val1_trawl1, val1_trawl2);
  
  arma::vec leb_mat0 = arma::trapz(supp0, val0);
  arma::vec leb_mat1 = arma::trapz(supp1, val1);
  
  double leb = leb_mat0(0) + leb_mat1(0);
  
  if (val1_trawl1(val1_trawl1.n_elem - 1) < val1_trawl2(val1_trawl2.n_elem - 1)) {
    leb += leb_AtA(arma::ones(1) * ub, trawl1, trawl_par1)(0);
  } else {
    leb += leb_AtA(arma::ones(1) * ub, trawl2, trawl_par2)(0);
  }
  
  return leb;
}
