#include "RcppArmadillo.h"
#include "levy_poisson.h"
#include "levy_skellam.h"
#include "levy_negbin.h"
#include "levy_dnegbin.h"
#include "levy_bivlog.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


double levy_intens(std::string levy_seed, arma::vec levy_par) {
  // wrapper function to compute Poisson intensity for some levy process
  //
  // arguments:
  // levy_seed  : levy seed of the process
  // levy_par   : corresponding parameters
  //
  // author: Dries Cornilly
  
  double intens = 0.0;
  if (levy_seed == "Poisson") {
    intens = intens_POISSON(levy_par(0));
  } else if (levy_seed == "Skellam") {
    intens = intens_SKELLAM(levy_par(0), levy_par(1));
  } else if (levy_seed == "negBin") {
    intens = intens_NEGBIN(levy_par(0), levy_par(1));
  } else if (levy_seed == "DnegBin") {
    intens = intens_DNEGBIN(levy_par(0), levy_par(1), levy_par(2), levy_par(3));
  } else if (levy_seed == "bivlog") {
    intens = intens_BIVLOG(levy_par(0), levy_par(1), levy_par(2));
  } else {
    stop("provide valid Lévy seed");
  }
  
  return intens;
}

arma::vec levy_rjump(int n, std::string levy_seed, arma::vec levy_par) {
  // wrapper function to compute jump sizes in a compound Poisson model
  //
  // arguments:
  // n          : sample size
  // levy_seed  : levy seed of the process
  // levy_par   : corresponding parameters
  //
  // author: Dries Cornilly
  
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

arma::mat levy_rjump_mv(int n, int k, std::string levy_seed, arma::vec levy_par) {
  // wrapper function to compute jump sizes in a multivariate compound Poisson model
  //
  // arguments:
  // n          : sample size
  // k          : dimension
  // levy_seed  : levy seed of the process
  // levy_par   : corresponding parameters
  //
  // author: Dries Cornilly
  
  arma::mat rj = arma::zeros(n, k);
  if (levy_seed == "bivlog") {
    double alpha1 = levy_par(1) / (1.0 - levy_par(1));
    double alpha2 = levy_par(2) / (1.0 - levy_par(2));
    rj = rjump_BIVLOG(n, alpha1 / (alpha1 + alpha2 + 1.0), alpha2 / (alpha1 + alpha2 + 1.0));
  } else {
    stop("provide valid Lévy seed");
  }
  
  return rj;
}

// [[Rcpp::export()]]
arma::vec levy_cum_fit(std::string levy_seed, double k1_sample, double k2_sample) {
  // wrapper function to fit univariate levy seed based on its first two cumulants
  //
  // arguments:
  // levy_seed  : levy seed of the process
  // k1_sample  : first cumulant
  // k2_sample  : second cumulant
  //
  // author: Dries Cornilly
  
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

// [[Rcpp::export()]]
arma::mat levy_varcovar(std::string levy_seed, arma::mat levy_par, arma::mat design_matrix) {
  // wrapper function to compute the variance covariance matrix of a multivariate levy seed,
  // see Veraart (2018)
  //
  // arguments:
  // levy_seed      : levy seed of the process
  // levy_par       : corresponding parameters
  // design_matrix  : design_matrix used in the multivariate levy seed
  //
  // author: Dries Cornilly
  
  arma::mat varcovar;
  if (levy_seed == "Poisson" || levy_seed == "Skellam") {
    varcovar = design_matrix * arma::diagmat(levy_par.col(0)) * design_matrix.t();
  } else if (levy_seed == "negBin") {
    int p = design_matrix.n_rows;
    int k = design_matrix.n_cols;
    
    varcovar = arma::zeros(p, p);
    for (int ii = 0; ii < k; ii++) {
      if (arma::sum(design_matrix.col(ii)) < 1.5) {
        double alpha = levy_par(ii, 1) / (1.0 - levy_par(ii, 1));
        double varV = levy_par(ii, 0) * alpha * alpha;
        double eV = levy_par(ii, 0) * alpha;
        
        arma::uvec ind = find(design_matrix.col(ii) > 0.5);
        varcovar(ind(0), ind(0)) += varV + eV;
      } else {
        double alpha_i = levy_par(ii, 1) / (1.0 - levy_par(ii, 1));
        double alpha_j = levy_par(ii, 2) / (1.0 - levy_par(ii, 2));
        
        double varU = levy_par(ii, 0);
        double eU = levy_par(ii, 0);
        
        arma::uvec ind = find(design_matrix.col(ii) > 0.5);
        varcovar(ind(0), ind(1)) += alpha_i * alpha_j * varU;
        varcovar(ind(1), ind(0)) += alpha_i * alpha_j * varU;
        varcovar(ind(0), ind(0)) += alpha_i * alpha_i * varU + alpha_i * eU;
        varcovar(ind(1), ind(1)) += alpha_j * alpha_j * varU + alpha_j * eU;
      } 
    }
  } else {
    stop("provide valid Lévy seed");
  }
  
  return varcovar;
}

