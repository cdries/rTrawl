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
