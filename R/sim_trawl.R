#' Simulation of Trawl processes
#'
#' simulates a path of a Trawl process
#'
#' TODO
#'
#' CITE TODO.
#' @name sim_trawl
#' @concept trawl
#' @param \dots any other passthru pareters
#' @author Dries Cornilly
#' @seealso \code{\link{fit_trawl}}
#' @references
#' TODO
#'
#' @examples
#'
#' TODO
#'
#' # simulations estimation
#' TODO
#'
#' @useDynLib rTrawl
#' @export sim_trawl
sim_trawl <- function(object, univariate = TRUE, ...) {
  
  # extract settings
  levy_seed <- object$levy_seed
  levy_par <- object$levy_par
  trawl <- object$trawl
  trawl_par <- object$trawl_par
  T0 <- object$T0
  TT <- object$TT
  observed_freq <- object$observed_freq
  b <- object$b
  
  # simulate
  if (univariate) {
    lsim <- simulate_trawl_uv(levy_seed, as.numeric(levy_par), trawl, as.numeric(trawl_par), 
                              as.numeric(T0), as.numeric(TT), as.numeric(observed_freq), as.numeric(b))
  } else {
    design_matrix <- object$design_matrix
    lsim <- simulate_trawl_mv(levy_seed, levy_par, trawl, trawl_par, design_matrix, as.numeric(T0), 
                              as.numeric(TT), as.numeric(observed_freq), as.numeric(b))
  }
  
  # export
  sim <- list("x_grid" = as.numeric(lsim$x_grid_observed),
              "p_grid" = as.numeric(lsim$p_grid_observed),
              "levy_seed" = levy_seed,
              "levy_par" = levy_par,
              "trawl" = trawl,
              "trawl_par" = trawl_par,
              "b" = b,
              "T0" = T0,
              "TT" = TT,
              "observed_freq" = observed_freq)
  
  if (hasArg(observe_latent)) {
    if (list(...)$observe_latent) {
      sim$x_grid <- as.numeric(lsim$x_grid_latent)
      sim$p_grid <- as.numeric(lsim$p_grid_latent)
    }
  }
  
  return(sim)
}