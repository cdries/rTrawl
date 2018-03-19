#' Estimation of Trawl processes
#'
#' estimates a path of a Trawl process
#'
#' TODO
#'
#' CITE TODO.
#' @name fit_levy
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
#' @export fit_levy
fit_levy <- function(object, ...) {
  
  # extract settings
  levy_seed <- object$levy_seed
  h <- object$h
  
  # trawl measure
  trawl <- object$trawl
  trawl_par <- object$trawl_par
  lebA <- leb_AtA(0.0, trawl, trawl_par)
  
  # cumulants
  x_grid <- object$x_grid
  p_grid <- object$p_grid
  T0 <- object$T0
  TT <- object$TT
  k1_sample <- cum_sample(1L, as.numeric(x_grid), as.numeric(p_grid), as.numeric(T0), as.numeric(TT)) / lebA
  k2_sample <- cum_sample(2L, as.numeric(x_grid), as.numeric(p_grid), as.numeric(T0), as.numeric(TT)) / lebA
  
  # fit using 1st (and 2nd) cumulant
  lfit <- levy_cum_fit(levy_seed, k1_sample, k2_sample)
  
  return(lfit)
}
