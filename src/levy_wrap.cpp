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

// [[Rcpp::export()]]
arma::vec levy_cum_fit(std::string levy_seed, double k1_sample, double k2_sample) {
  
  arma::vec levy_par;
  if (levy_seed == "Poisson") {
    levy_par = fit_POISSON(k1_sample);
  } else if (levy_seed == "Skellam") {
    levy_par = fit_SKELLAM(k1_sample, k2_sample);
  } else if (levy_seed == "negBin") {
    levy_par = fit_NEGBIN(k1_sample, k2_sample);
  // } else if (levy_seed == "DnegBin") {
  //   levy_par = fit_DNEGBIN(k1_sample, k2_sample);
  } else {
    stop("provide valid Lévy seed");
  }
  
  return levy_par;
}

arma::mat levy_varcovar(std::string levy_seed, arma::mat levy_par, arma::mat design_matrix) {
  
  arma::mat varcovar;
  if (levy_seed == "Poisson" || levy_seed == "Skellam") {
    varcovar = design_matrix * arma::diagmat(levy_par.col(0)) * design_matrix.t();
  // } else if (levy_seed == "negBin") {
  //   // TODO
  } else {
    stop("prove valid Lévy seed");
  }
  
  return varcovar;
}

