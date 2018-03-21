#include "RcppArmadillo.h"
#include "trawl_exp.h"
#include "trawl_gamma.h"
#include "trawl_invgauss.h"
#include "trawl_gig.h"
#include "trawl_helper.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


arma::vec trawl_times(arma::vec unif_seed, std::string trawl, arma::vec trawl_par,
                      double observed_freq, double Tmax, double b) {
  
  arma::vec surv_times(unif_seed.n_elem);
  if (trawl == "exp") {
    surv_times = survival_EXP(unif_seed, trawl_par(0), Tmax, b);
  } else if (trawl == "gamma") {
    surv_times = survival_GAMMA(unif_seed, trawl_par(0), trawl_par(1), Tmax, b);
  } else if (trawl == "invGauss") {
    surv_times = survival_INVGAUSS(unif_seed, trawl_par(0), trawl_par(1), Tmax, b, observed_freq);
  } else if (trawl == "gig") {
    surv_times = survival_GIG(unif_seed, trawl_par(0), trawl_par(1), trawl_par(2), Tmax, b, observed_freq);
  } else {
    stop("provide a valid trawl");
  }
  
  return surv_times;
}

// [[Rcpp::export()]]
int number_parameters_trawl(std::string trawl) {
  
  int n = 0;
  if (trawl == "exp") {
    n = 1;
  } else if (trawl == "gamma") {
    n = 2;
  } else if (trawl == "invGauss") {
    n = 2;
  } else if (trawl == "gig") {
    n = 3;
  } else {
    stop("provide a valid trawl");
  }
  
  return n;
}

// [[Rcpp::export()]]
List trawl_bounds(std::string trawl) {
  
  arma::vec lb;
  arma::vec ub;
  if (trawl == "exp") {
    lb = ub = arma::ones(1) * 1e-7;
    ub(0) = arma::datum::inf;
  } else if (trawl == "gamma") {
    lb = ub = arma::ones(2) * 1e-7;
    lb(1) += 1.0;
    ub(0) = ub(1) = arma::datum::inf;
  } else if (trawl == "invGauss") {
    lb = ub = arma::ones(2) * 1e-7;
    ub(0) = ub(1) = arma::datum::inf;
  } else if (trawl == "gig") {
    lb = ub = arma::ones(3) * 1e-7;
    lb(2) = -arma::datum::inf;
    ub(0) = ub(1) = ub(2) = arma::datum::inf;
  } else {
    stop("provide a valid trawl");
  }
  
  List bounds;
  bounds["lb"] = lb;
  bounds["ub"] = ub;
  
  return bounds;
}

// [[Rcpp::export()]]
arma::vec trawl_x0(std::string trawl) {
  
  arma::vec x0;
  if (trawl == "exp") {
    x0 = arma::ones(1) * 2.0;
  } else if (trawl == "gamma") {
    x0 = arma::ones(2) * 0.2;
    x0(1) = 1.9;
  } else if (trawl == "invGauss") {
    x0 = arma::ones(2) * 0.1;
    x0(1) = 0.5;
  } else if (trawl == "gig") {
    x0 = arma::ones(3) * 0.1;
    x0(1) = 0.9;
    x0(2) = -0.2;
  } else {
    stop("provide a valid trawl");
  }
  
  return x0;
}

// [[Rcpp::export()]]
arma::vec leb_AtA(arma::vec h, std::string trawl, arma::vec trawl_par) {
  
  arma::vec leb;
  if (trawl == "exp") {
    leb = leb_AtA_EXP(h, trawl_par(0));
  } else if (trawl == "gamma") {
    leb = leb_AtA_GAMMA(h, trawl_par(0), trawl_par(1));
  } else if (trawl == "invGauss") {
    leb = leb_AtA_INVGAUSS(h, trawl_par(0), trawl_par(1));
  } else if (trawl == "gig") {
    leb = leb_AtA_GIG(h, trawl_par(0), trawl_par(1), trawl_par(2));
  } else {
    stop("provide a valid trawl");
  }
  
  return leb;
}

arma::mat d_leb_AtA(arma::vec h, std::string trawl, arma::vec trawl_par) {
  
  arma::mat d_leb;
  if (trawl == "exp") {
    d_leb = d_leb_AtA_EXP(h, trawl_par(0));
  } else if (trawl == "gamma") {
    d_leb = d_leb_AtA_GAMMA(h, trawl_par(0), trawl_par(1));
  } else if (trawl == "invGauss") {
    d_leb = d_leb_AtA_INVGAUSS(h, trawl_par(0), trawl_par(1));
  } else if (trawl == "gig") {
    d_leb = d_leb_AtA_GIG(h, trawl_par(0), trawl_par(1), trawl_par(2));
  } else {
    stop("provide a valid trawl");
  }
  
  return d_leb;
}

arma::vec leb_autocorrelator(arma::vec h, std::string trawl1, arma::vec trawl_par1,
                             std::string trawl2, arma::vec trawl_par2) {
  
  int n_h = h.n_elem;
  arma::vec leb = arma::zeros(n_h);
  for (int ii = 0; ii < n_h; ii++) {
    
    double hii = h(ii);
    std::string trawl1_ii = trawl1;
    std::string trawl2_ii = trawl2;
    arma::vec trawl_par1_ii = trawl_par1;
    arma::vec trawl_par2_ii = trawl_par2;
    if (hii < 0) {
      hii = -hii;
      trawl1_ii = trawl2;
      trawl2_ii = trawl1;
      trawl_par1_ii = trawl_par2;
      trawl_par2_ii = trawl_par1;
    }
    
    if (trawl1 == "exp" && trawl2 == "exp") {
      leb(ii) = leb_EXP_EXP(hii, trawl_par1_ii(0), trawl_par2_ii(0));
    } else {
      leb(ii) = leb_GEN_GEN(hii, trawl1_ii, trawl_par1_ii, trawl2_ii, trawl_par2_ii);
    }
  }
  
  return leb;
}

// [[Rcpp::export()]]
arma::vec trawl_function(arma::vec h, std::string trawl, arma::vec trawl_par) {
  
  arma::vec val;
  if (trawl == "exp") {
    val = trawl_EXP(h, trawl_par(0));
  } else if (trawl == "gamma") {
    val = trawl_GAMMA(h, trawl_par(0), trawl_par(1));
  } else if (trawl == "invGauss") {
    val = trawl_INVGAUSS(h, trawl_par(0), trawl_par(1));
  } else if (trawl == "gig") {
    val = trawl_GIG(h, trawl_par(0), trawl_par(1), trawl_par(2));
  } else {
    stop("provide a valid trawl");
  }
  
  return val;
}

