#' Simulation of Trawl processes
#'
#' simulates a path of a Trawl process
#'
#' TODO
#'
#' CITE TODO.
#' @name sim_trawl
#' @concept trawl
#' @param object bla
#' @param univariate bla
#' @param \dots any other passthru pareters
#' @author Dries Cornilly
#' @seealso \code{\link{fit_trawl}}
#' @references
#' TODO
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
