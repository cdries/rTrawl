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
  // compute the time an observation is coverd by the trawl set 
  //
  // arguments:
  // unif_seed  : standard uniform randomly generated numbers
  // trawl      : trawl of the process
  // trawl_par  : trawl parameters
  // observed_freq : determines how accurate survival times should be approximated
  // Tmax       : maximum observation window, times larger than this are not useful
  // b          : b parameters, as in Shephard and Yang (2017)
  //
  // author: Dries Cornilly
  
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
  // gives the number of parameters for each trawl
  //
  // arguments:
  // trawl      : trawl of the process
  //
  // author: Dries Cornilly
  
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
  // provides upper and lower bounds for each of the trawl parameters
  //
  // arguments:
  // trawl      : trawl of the process
  //
  // author: Dries Cornilly
  
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
  // provides initial values for each of the trawl parameters
  //
  // arguments:
  // trawl      : trawl of the process
  //
  // author: Dries Cornilly
  
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
  // wrapper function to compute the lebesgue measure of A_0 cap A_h
  //
  // arguments:
  // h          : determines A_h, if 0, then leb(A) is computed
  // trawl      : trawl of the process
  // trawl_par  : trawl parameters
  //
  // author: Dries Cornilly
  
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
  // wrapper function to compute gradient the lebesgue measure of A_0 cap A_h
  // with respect to the trawl parameters
  //
  // arguments:
  // h          : determines A_h, if 0, then leb(A) is computed
  // trawl      : trawl of the process
  // trawl_par  : trawl parameters
  //
  // author: Dries Cornilly
  
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
  // wrapper function to compute the lebesgue measure of the minimum of trawl 1 and trawl 2
  //
  // arguments:
  // h          : argument for A_h^1 cap A_0^2 with h positive and -h at 2 for h negative
  // trawl1     : trawl 1
  // trawl_par1 : trawl parameters for trawl 1
  // trawl2     : trawl 2
  // trawl_par2 : trawl parameters for trawl 2
  //
  // author: Dries Cornilly
  
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
  // wrapper function to each of the trawl functions
  //
  // arguments:
  // h          : argument at which to compute the function
  // trawl1     : trawl 1
  // trawl_par1 : trawl parameters for trawl 1
  // trawl2     : trawl 2
  // trawl_par2 : trawl parameters for trawl 2
  //
  // author: Dries Cornilly
  
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

double leb_autocorrelator_general(double t1, double t2, double s1, double s2, 
                                  double b1, double b2, bool area1, bool area2,
                                  std::string trawl1, std::string trawl2,
                                  arma::vec trawl_par1, arma::vec trawl_par2) {
  // helper function to compute the lebesgue measure at the intersection of two 
  // trawls, mainly used when working with differenced series
  //
  // arguments:
  // t1         : lower bound for trawl 1
  // t2         : upper bound for trawl 1
  // s1         : lower bound for trawl 2
  // s2         : upper bound for trawl 2
  // b1         : b parameter for trawl 1, see Shephard and Yang (2017)
  // b2         : b parameter for trawl 2, see Shephard and Yang (2017)
  // area1      : boolean indicating if the area under b1 is measured or not for trawl 1
  // area2      : boolean indicating if the area under b2 is measured or not for trawl 2
  // trawl1     : trawl 1
  // trawl2     : trawl 2
  // trawl_par1 : trawl parameters for trawl 1
  // trawl_par2 : trawl parameters for trawl 2
  //
  // author: Dries Cornilly
  
  double leb = 0.0;
  double Nsupp1 = 100000;
  
  double TT = s2;
  if (t2 < s2) TT = t2;
  
  if (area1 && area2) {
    double T0 = s1;
    if (t1 > s1) T0 = t1;
    
    if (T0 < TT) {
      arma::vec supp = -TT + (arma::logspace(0.0, 1.0, Nsupp1) - 1.0) * (TT - T0) / 9.0;
      
      arma::vec val_trawl1 = b1 + (1.0 - b1) * trawl_function(-supp - t2, trawl1, trawl_par1);
      arma::vec val_trawl2 = b2 + (1.0 - b2) * trawl_function(-supp - s2, trawl2, trawl_par2);
      
      arma::vec val = arma::min(val_trawl1, val_trawl2);
      arma::mat tmp = arma::trapz(supp, val);
      leb = tmp(0);
    }
  } else {
    double T0 = TT - 100.0;
    if (area1) T0 = t1;
    if (area2) T0 = s1;
    
    if (T0 < TT) {
      arma::vec supp = -TT + (arma::logspace(0.0, 1.0, Nsupp1) - 1.0) * (TT - T0) / 9.0;
      
      arma::vec val_trawl1 = b1 + (1.0 - b1) * trawl_function(-supp - t2, trawl1, trawl_par1);
      arma::vec val_trawl2 = b2 + (1.0 - b2) * trawl_function(-supp - s2, trawl2, trawl_par2);
      
      arma::vec val_trawl1_min = arma::zeros(Nsupp1);
      arma::vec val_trawl2_min = arma::zeros(Nsupp1);
      if (!area1) {
        val_trawl1 = b1 + (1.0 - b1) * trawl_function(-supp - t1, trawl1, trawl_par1);
        val_trawl1_min = b1 + (1.0 - b1) * trawl_function(-supp - t2, trawl1, trawl_par1);
      }
      if (!area2) {
        val_trawl2 = b2 + (1.0 - b2) * trawl_function(-supp - s1, trawl2, trawl_par2);
        val_trawl2_min = b2 + (1.0 - b2) * trawl_function(-supp - s2, trawl2, trawl_par2);
      }
      arma::vec val_min = arma::max(val_trawl1_min, val_trawl2_min);
      
      arma::vec val = arma::min(val_trawl1, val_trawl2);
      arma::uvec ind = find(val < val_min);
      val.elem(ind) = val_min.elem(ind);
      arma::mat tmp = arma::trapz(supp, val - val_min);
      
      leb = tmp(0);
    }
  }
  
  return leb;
}
