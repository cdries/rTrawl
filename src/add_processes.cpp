#include "RcppArmadillo.h"
#include "observe_process.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export()]]
List add_processes(arma::vec x_grid1, arma::vec p_grid1, arma::vec x_grid2, arma::vec p_grid2,
                   double T0, double TT, double observed_freq) {
  // returns the sum of two non-equidistant time series at a chosen observation frequency,
  // handles non-equidistant observations and does not return an equidistant grid
  //
  // arguments:
  // x_grid1        : vector with observed times for process 1
  // p_grid1        : corresponding process values for process 1
  // x_grid2        : vector with observed times for process 2
  // p_grid2        : corresponding process values for process 2
  // T0             : beginpoint of observation interval
  // TT             : endpoint of observation interval
  // observed_freq  : new observation frequency
  //
  // author: Dries Cornilly
  
  // initialize
  int n1 = x_grid1.n_elem;
  int n2 = x_grid2.n_elem;
  int ii1 = 0;
  int ii2 = 0;
  arma::vec diff1 = arma::diff(p_grid1);
  arma::vec diff2 = arma::diff(p_grid2);
  
  // find starting point of sum process
  if (x_grid1(0) < x_grid2(0)) {
    while (x_grid1(ii1) < x_grid2(0)) ii1++;
  } else {
    while (x_grid1(0) > x_grid2(ii2)) ii2++;
  }
  
  // initialize sum process
  arma::vec x_grid_joint = arma::zeros(n1 + n2 - ii1 - ii2);
  arma::vec p_grid_joint = arma::zeros(n1 + n2 - ii1 - ii2);
  int n_joint = x_grid_joint.n_elem;
  
  // initial value of the process
  double current_value;
  if (ii1 > ii2) {
    
    // grid 1 jumped over grid 2 when finding starting point
    current_value = p_grid1(ii1 - 1) + p_grid2(ii2);
    x_grid_joint(0) = x_grid2(ii2);
    ii2++;
    
  } else {
    
    // grid 2 jumped over grid 1 when finding starting point
    current_value = p_grid2(ii2 - 1) + p_grid1(ii1);
    x_grid_joint(0) = x_grid1(ii1);
    ii1++;
  }
  p_grid_joint(0) = current_value;
  
  // iterate over the time grids
  int jj = 1;
  while (ii1 < n1 && ii2 < n2) {

    if (x_grid1(ii1) < x_grid2(ii2)) {

      // fill in from series 1
      x_grid_joint(jj) = x_grid1(ii1);
      current_value += diff1(ii1 - 1);
      p_grid_joint(jj) = current_value;
      ii1++;

    } else {

      // fill in from series 2
      x_grid_joint(jj) = x_grid2(ii2);
      current_value += diff2(ii2 - 1);
      p_grid_joint(jj) = current_value;
      ii2++;

    }

    jj++;

  }
  
  // finish up last part
  if (ii1 >= n1) {

    // add last part of series 2
    x_grid_joint.subvec(jj, n_joint - 1) = x_grid2.subvec(ii2, n2 - 1);
    p_grid_joint.subvec(jj, n_joint - 1) = current_value + diff2.subvec(ii2 - 1, n2 - 2);

  } else {

    // add last part of series 1
    x_grid_joint.subvec(jj, n_joint - 1) = x_grid1.subvec(ii1, n1 - 1);
    p_grid_joint.subvec(jj, n_joint - 1) = current_value + diff1.subvec(ii1 - 1, n1 - 2);
  }

  // return list  
  List out = observe_process(x_grid_joint, p_grid_joint, T0, TT, observed_freq);
  StringVector outnames(2);
  outnames(0) = "x_grid";
  outnames(1) = "p_grid";
  out.names() = outnames;
  out["T0"] = T0;
  out["TT"] = TT;
  
  return out;
}
