#include "RcppArmadillo.h"
#include "levy_wrap.h"
#include "trawl_wrap.h"
#include "observe_process.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export()]]
List simulate_trawl_uv(std::string levy_seed, arma::vec levy_par, std::string trawl, 
                       arma::vec trawl_par, double T0, double TT, double observed_freq, double b) {
  
  // latent time grid
  double intens = levy_intens(levy_seed, levy_par);
  int n_jumps = rpois(1, 1.5 * (TT - T0) * intens)(0);
  arma::vec jump_times = runif(n_jumps, T0 - 0.5 * (TT - T0), TT);
  arma::vec x_grid_latent = arma::sort(jump_times);
  
  // survival times and jump sizes
  arma::vec unif_seed = runif(n_jumps, 0.0, 1.0);
  arma::vec survival_times = trawl_times(unif_seed, trawl, trawl_par, observed_freq, 2 * (TT - T0), b);
  arma::vec jump_sizes = levy_rjump(n_jumps, levy_seed, levy_par);
  
  // latent trawl process
  x_grid_latent = arma::join_cols(x_grid_latent, x_grid_latent + survival_times);
  arma::vec p_grid_latent = arma::join_cols(jump_sizes, -jump_sizes);
  arma::uvec order_x_grid = arma::sort_index(x_grid_latent);
  x_grid_latent = x_grid_latent.elem(order_x_grid);
  p_grid_latent = arma::cumsum(p_grid_latent.elem(order_x_grid));
  arma::uvec ind = find((x_grid_latent >= T0) && (x_grid_latent < TT));
  x_grid_latent = x_grid_latent.elem(ind);
  p_grid_latent = p_grid_latent.elem(ind);
  
  // observe process
  List process_observed = observe_process(x_grid_latent, p_grid_latent, T0, TT, observed_freq);
  
  // return list
  List out;
  out["x_grid_latent"] = x_grid_latent;
  out["p_grid_latent"] = p_grid_latent;
  out["x_grid_observed"] = process_observed["x_grid_observed"];
  out["p_grid_observed"] = process_observed["p_grid_observed"];
  
  return out;
}

// [[Rcpp::export()]]
List simulate_trawl_mv_Poisson(arma::mat levy_par, List trawl, List trawl_par, 
                               arma::mat design_matrix, double T0, double TT, 
                               double observed_freq, arma::vec b) {
  
  int p = design_matrix.n_rows;
  int k = design_matrix.n_cols;
  List x_grid_latent(p);
  List jump_sizes(p);
  List survival_times(p);
  arma::vec components = arma::zeros(p);
  
  // iterate over the latent factors
  for (int ii = 0; ii < k; ii++) {
    
    // latent time grid
    double intens = levy_intens("Poisson", levy_par.row(ii));
    int n_jumps = rpois(1, 1.5 * (TT - T0) * intens)(0);
    arma::vec jump_times = runif(n_jumps, T0 - 0.5 * (TT - T0), TT);
    arma::vec x_grid_latent_factor = arma::sort(jump_times);
    
    // jump sizes and random seed
    arma::vec jump_sizes_factor = levy_rjump(n_jumps, "Poisson", levy_par.row(ii));
    arma::vec unif_seed = runif(n_jumps, 0.0, 1.0);
    
    // combine and add jump sizes
    for (int jj = 0; jj < p; jj++) {
      
      if (design_matrix(jj, ii) * design_matrix(jj, ii) > 0.5) {
        
        if (components(jj) < 0.5) { // if list is still empty
          components(jj)++;
          
          x_grid_latent(jj) = x_grid_latent_factor;
          
          survival_times(jj) = trawl_times(unif_seed, trawl(jj), trawl_par(jj),
                         observed_freq, 2 * (TT - T0), b(jj));
          
          if (design_matrix(jj, ii) > 0.5) {
            jump_sizes(jj) = jump_sizes_factor;
          } else {
            jump_sizes(jj) = -jump_sizes_factor;
          }
        } else { // if we want to append to current lists
          arma::vec x_grid_latent_temp = x_grid_latent(jj);
          x_grid_latent(jj) = arma::join_cols(x_grid_latent_temp, x_grid_latent_factor);
          
          arma::vec survival_times_temp = survival_times(jj);
          survival_times(jj) = arma::join_cols(survival_times_temp,
                         trawl_times(unif_seed, trawl(jj), trawl_par(jj),
                                     observed_freq, 2 * (TT - T0), b(jj)));
          
          arma::vec jump_sizes_temp = jump_sizes(jj);
          if (design_matrix(jj, ii) > 0.5) {
            jump_sizes(jj) = arma::join_cols(jump_sizes_temp, jump_sizes_factor);
          } else {
            jump_sizes(jj) = arma::join_cols(jump_sizes_temp, -jump_sizes_factor);
          }
        }
      }
    }
  }
  
  // latent trawl process - observe process
  List p_grid_latent(p);
  List x_grid_observed(p);
  List p_grid_observed(p);
  for (int ii = 0; ii < p; ii++) {
    
    // latent
    arma::vec x_grid_temp = x_grid_latent(ii);
    arma::vec survival_times_temp = survival_times(ii);
    arma::vec jump_sizes_temp = jump_sizes(ii);
    
    x_grid_temp = arma::join_cols(x_grid_temp, x_grid_temp + survival_times_temp);
    arma::vec p_grid_temp = arma::join_cols(jump_sizes_temp, -jump_sizes_temp);
    arma::uvec order_x_grid = arma::sort_index(x_grid_temp);
    x_grid_temp = x_grid_temp.elem(order_x_grid);
    p_grid_temp = arma::cumsum(p_grid_temp.elem(order_x_grid));
    arma::uvec ind = find((x_grid_temp >= T0) && (x_grid_temp < TT));
    x_grid_temp = x_grid_temp.elem(ind);
    p_grid_temp = p_grid_temp.elem(ind);
    
    x_grid_latent(ii) = x_grid_temp;
    p_grid_latent(ii) = p_grid_temp;
    
    // observe
    List process_observed = observe_process(x_grid_temp, p_grid_temp, T0, TT, observed_freq);
    x_grid_observed(ii) = process_observed["x_grid_observed"];
    p_grid_observed(ii) = process_observed["p_grid_observed"];
  }
  
  // return list
  List out;
  out["x_grid_latent"] = x_grid_latent;
  out["p_grid_latent"] = p_grid_latent;
  out["x_grid_observed"] = x_grid_observed;
  out["p_grid_observed"] = p_grid_observed;
  
  return out;
}

// [[Rcpp::export()]]
List simulate_trawl_mv_negBin(arma::mat levy_par, List trawl, List trawl_par, 
                              arma::mat design_matrix, double T0, double TT, 
                              double observed_freq, arma::vec b) {
  
  int p = design_matrix.n_rows;
  int k = design_matrix.n_cols;
  List x_grid_latent(p);
  List jump_sizes(p);
  List survival_times(p);
  arma::vec components = arma::zeros(p);
  
  // iterate over the latent factors
  for (int ii = 0; ii < k; ii++) {
    
    // latent time grid
    double intens = 0.0;
    if (arma::sum(design_matrix.col(ii)) < 1.5) {
      intens = levy_intens("negBin", levy_par.row(ii).t());
    } else {
      intens = levy_intens("bivlog", levy_par.row(ii).t());
    }
    int n_jumps = rpois(1, 1.5 * (TT - T0) * intens)(0);
    arma::vec jump_times = runif(n_jumps, T0 - 0.5 * (TT - T0), TT);
    arma::vec x_grid_latent_factor = arma::sort(jump_times);
    
    // jump sizes and random seed
    arma::mat jump_sizes_factor;
    if (arma::sum(design_matrix.col(ii)) < 1.5) {
      jump_sizes_factor = levy_rjump(n_jumps, "negBin", levy_par.row(ii).t());
    } else {
      jump_sizes_factor = levy_rjump_mv(n_jumps, 2, "bivlog", levy_par.row(ii).t());
    }
    arma::vec unif_seed = runif(n_jumps, 0.0, 1.0);
    
    // combine and add jump sizes
    int kk = 0;
    for (int jj = 0; jj < p; jj++) {
      
      if (design_matrix(jj, ii) > 0.5) {
        
        if (components(jj) < 0.5) { // if list is still empty
          components(jj)++;
          
          x_grid_latent(jj) = x_grid_latent_factor;
          
          survival_times(jj) = trawl_times(unif_seed, trawl(jj), trawl_par(jj),
                         observed_freq, 2 * (TT - T0), b(jj));
          
          jump_sizes(jj) = jump_sizes_factor.col(kk);
          kk++;
        } else { // if we want to append to current lists
          arma::vec x_grid_latent_temp = x_grid_latent(jj);
          x_grid_latent(jj) = arma::join_cols(x_grid_latent_temp, x_grid_latent_factor);
          
          arma::vec survival_times_temp = survival_times(jj);
          survival_times(jj) = arma::join_cols(survival_times_temp,
                         trawl_times(unif_seed, trawl(jj), trawl_par(jj),
                                     observed_freq, 2 * (TT - T0), b(jj)));
          
          arma::vec jump_sizes_temp = jump_sizes(jj);
          jump_sizes(jj) = arma::join_cols(jump_sizes_temp, jump_sizes_factor.col(kk));
          kk++;
        }
      }
    }
  }
  
  // latent trawl process - observe process
  List p_grid_latent(p);
  List x_grid_observed(p);
  List p_grid_observed(p);
  for (int ii = 0; ii < p; ii++) {
    
    // latent
    arma::vec x_grid_temp = x_grid_latent(ii);
    arma::vec survival_times_temp = survival_times(ii);
    arma::vec jump_sizes_temp = jump_sizes(ii);
    
    x_grid_temp = arma::join_cols(x_grid_temp, x_grid_temp + survival_times_temp);
    arma::vec p_grid_temp = arma::join_cols(jump_sizes_temp, -jump_sizes_temp);
    arma::uvec order_x_grid = arma::sort_index(x_grid_temp);
    x_grid_temp = x_grid_temp.elem(order_x_grid);
    p_grid_temp = arma::cumsum(p_grid_temp.elem(order_x_grid));
    arma::uvec ind = find((x_grid_temp >= T0) && (x_grid_temp < TT));
    x_grid_temp = x_grid_temp.elem(ind);
    p_grid_temp = p_grid_temp.elem(ind);
    
    x_grid_latent(ii) = x_grid_temp;
    p_grid_latent(ii) = p_grid_temp;
    
    // observe
    List process_observed = observe_process(x_grid_temp, p_grid_temp, T0, TT, observed_freq);
    x_grid_observed(ii) = process_observed["x_grid_observed"];
    p_grid_observed(ii) = process_observed["p_grid_observed"];
  }
  
  // return list
  List out;
  out["x_grid_latent"] = x_grid_latent;
  out["p_grid_latent"] = p_grid_latent;
  out["x_grid_observed"] = x_grid_observed;
  out["p_grid_observed"] = p_grid_observed;
  
  return out;
}
