#' Simulation of trawl processes
#'
#' simulation of univariate or multivariate trawl processes
#'
#' the passed object contains all the process specifications in list format. The following
#' arguments should ideally be present in the object: levy seed (levy_seed), corresponding
#' levy parameters (levy_par), trawl (trawl), corresponding trawl parameters (trawl_par),
#' initial observation time (T0), end of observation period (TT), the frequency at which
#' the process is observed (observed_freq), and the b parameter governing the distinction
#' between levy and trawl process changes. In case of multivariate simulation, also the
#' design matrix showing the underlying factors should be provided (design_matrix).
#' 
#' The observation frequency is for example 1 second, 
#' which means that a change in the process value is only observed the next round second. 
#' For the b parameter, we refer to Shephard and Yang (2017). Some default parameters are set,
#' but it is best not to rely on these. See examples.
#'
#' @name sim_trawl
#' @concept trawl
#' @param object object containing all the specifications for the process, see details
#' @param univariate boolean indicating univariate or multivariate simulation
#' @param \dots any other passthrough parameters
#' @return x_grid: vector of time points where the process value is given
#' @return p_grid: vector with process values for the time points in x_grid
#' @return the other return arguments echo the process specifications
#' @author Dries Cornilly
#' @seealso \code{\link{fit_trawl}}
#' @references
#' Barndorff‐Nielsen, O. E., Lunde, A., Shephard, N., & Veraart, A. E. (2014). 
#' Integer‐valued Trawl Processes: A Class of Stationary Infinitely Divisible Processes. 
#' Scandinavian Journal of Statistics, 41(3), 693-724.
#' 
#' Shephard, N., & Yang, J. J. (2017). 
#' Continuous time analysis of fleeting discrete price moves. 
#' Journal of the American Statistical Association, 112(519), 1090-1106.
#' 
#' Veraart, A. E. (2018). 
#' Modelling, simulation and inference for multivariate time series of counts. 
#' arXiv preprint arXiv:1608.03154.
#'
#' @examples
#'
#' #TODO
#'
#' # simulations estimation
#' #TODO
#'
#' @importFrom methods hasArg
#' @importFrom Rcpp evalCpp
#' @useDynLib rTrawl
#' @export sim_trawl
sim_trawl <- function(object, univariate = TRUE, ...) {
  
  # extract settings
  nn <- names(object)
  if ("levy_seed" %in% nn) levy_seed <- object$levy_seed else levy_seed <- "Poisson"
  if ("levy_par" %in% nn) levy_par <- object$levy_par else levy_par <- 0.0131
  if ("trawl" %in% nn) trawl <- object$trawl else trawl <- "exp"
  if ("trawl" %in% nn) trawl_par <- object$trawl_par else trawl_par <- trawl_x0(trawl)
  if ("T0" %in% nn) T0 <- object$T0 else T0 <- 0
  if ("TT" %in% n) TT <- object$TT else TT <- 3600
  if ("observed_freq" %in% nn) observed_freq <- object$observed_freq else observed_freq <- 1e-3
  if ("b" %in% n) b <- object$b else b <- 0
  
  
  # simulate
  if (univariate) {
    lsim <- simulate_trawl_uv(levy_seed, as.numeric(levy_par), trawl, as.numeric(trawl_par), 
                              as.numeric(T0), as.numeric(TT), as.numeric(observed_freq), as.numeric(b))
  } else {
    if ("design_matrix" %in% nn) design_matrix <- object$design_matrix else design_matrix <- diag(length(trawl))
    
    if (levy_seed %in% c("Poisson", "Skellam")) {
      lsim <- simulate_trawl_mv_Poisson(levy_par, trawl, trawl_par, design_matrix, as.numeric(T0), 
                                        as.numeric(TT), as.numeric(observed_freq), as.numeric(b))
    } else if (levy_seed == "negBin") {
      lsim <- simulate_trawl_mv_negBin(levy_par, trawl, trawl_par, design_matrix, as.numeric(T0), 
                                       as.numeric(TT), as.numeric(observed_freq), as.numeric(b))
    } else {
      stop("Multivariate trawls are not implemented for this choice of Levy seed.")
    }
  }
  
  # export
  sim <- list("x_grid" = lsim$x_grid_observed,
              "p_grid" = lsim$p_grid_observed,
              "levy_seed" = levy_seed,
              "levy_par" = levy_par,
              "trawl" = trawl,
              "trawl_par" = trawl_par,
              "b" = b,
              "T0" = T0,
              "TT" = TT,
              "observed_freq" = observed_freq)
  
  if (hasArg(observe_latent)) {
    observe_latent <- list(...)$observe_latent
    if (observe_latent) {
      sim$x_grid <- lsim$x_grid_latent
      sim$p_grid <- lsim$p_grid_latent
    }
  }
  
  return(sim)
}
