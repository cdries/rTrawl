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
  
  // initialize joint process
  int n1 = x_grid1.n_elem;
  int n2 = x_grid2.n_elem;
  arma::vec x_grid_joint = arma::zeros(n1 + n2);
  arma::vec p_grid_joint = arma::zeros(n1 + n2);
  
  // iterate over the time grids
  int ii1 = 0;
  int ii2 = 0;
  double current_value = 0.0;
  while (ii1 < n1 && ii2 < n2) {
    
    if (x_grid1(ii1) < x_grid2(ii2)) {
      
      // fill in from series 1
      current_value += p_grid1(ii1);
      x_grid_joint(ii1 + ii2) = current_value;
      ii1++;
      
    } else {
      
      // fill in from series 2
      current_value += p_grid2(ii2);
      x_grid_joint(ii1 + ii2) = current_value;
      ii2++;
      
    }
  }
  
  // finish up last part
  if (ii1 >= n1) {
    
    // add last part of series 2
    x_grid_joint.subvec(ii1 + ii2, n1 + n2 - 2) = current_value + p_grid2.subvec(ii2, n2 - 1);
    
  } else {
    
    // add last part of series 1
    x_grid_joint.subvec(ii1 + ii2, n1 + n2 - 2) = current_value + p_grid1.subvec(ii1, n1 - 1);
  }
  
  // return list
  List out = observe_process(x_grid_joint, p_grid_joint, T0, TT, observed_freq);
  
  return out;
}
