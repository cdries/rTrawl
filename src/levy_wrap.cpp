#include "RcppArmadillo.h"
#include "levy_poisson.h"
#include "levy_skellam.h"
#include "levy_negbin.h"
#include "levy_dnegbin.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


double levy_intens(std::string levy_seed, arma::vec levy_par) {
  
  double intens = 0.0;
  if (levy_seed == "Poisson") {
    intens = intens_POISSON(levy_par(0));
  } else if (levy_seed == "Skellam") {
    intens = intens_SKELLAM(levy_par(0), levy_par(1));
  } else if (levy_seed == "negBin") {
    intens = intens_NEGBIN(levy_par(0), levy_par(1));
  } else if (levy_seed == "DnegBin") {
    intens = intens_DNEGBIN(levy_par(0), levy_par(1), levy_par(2), levy_par(3));
  } else {
    stop("provide valid Lévy seed");
  }
  
  return intens;
}

arma::vec levy_rjump(int n, std::string levy_seed, arma::vec levy_par) {
  
  arma::vec rj = arma::zeros(n);
  if (levy_seed == "Poisson") {
    rj = rjump_POISSON(n);
  } else if (levy_seed == "Skellam") {
    rj = rjump_SKELLAM(n, levy_par(0), levy_par(1));
  } else if (levy_seed == "negBin") {
    rj = rjump_NEGBIN(n, levy_par(1));
  } else if (levy_seed == "DnegBin") {
    rj = rjump_DNEGBIN(n, levy_par(0), levy_par(1), levy_par(2), levy_par(3));
  } else {
    stop("provide valid Lévy seed");
  }
  
  return rj;
}
