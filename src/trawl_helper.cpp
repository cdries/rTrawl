#include "RcppArmadillo.h"
#include "trawl_exp.h"
#include "trawl_gamma.h"
#include "trawl_invgauss.h"
#include "trawl_gig.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// h always positive!
double leb_EXP_EXP(double h, double lambda1, double lambda2) {
  
  double leb;
  
  if (lambda2 < lambda1) {
    leb = exp(-lambda2 * h) / lambda2 - 
      (1 / lambda2 - 1 / lambda1) * exp(lambda1 * lambda2 * h / (lambda2 - lambda1));
  } else {
    leb = exp(-lambda2 * h) / lambda2;
  }
  
  return leb; 
}

// h always positive!
double leb_GEN_GEN(double  h, std::string trawl1, arma::vec trawl_par1, 
                   std::string trawl2, arma::vec trawl_par2) {
  
  double leb = 0.0;
  
  return leb;
}
