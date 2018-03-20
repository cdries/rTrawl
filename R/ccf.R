#' Sample ACF of irregularly observed time series
#'
#' sample autocorrelation of a trawl process
#'
#' TODO
#'
#' CITE TODO.
#' @name ccf_sample
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
#' @export ccf_sample
#' @useDynLib rTrawl
ccf_sample <- function(object, h, dff = 0, lag_max = 25, ...) {
  
  x_grid1 <- object$x_grid[[1]]
  x_grid2 <- object$x_grid[[2]]
  p_grid1 <- object$p_grid[[1]]
  p_grid2 <- object$p_grid[[2]]
  TT <- object$TT
  
  if (dff < 0.5) {
    # cross correlation function of the process itself
    ccfh <- as.numeric(ccf_sample_p(h, x_grid1, p_grid1, x_grid2, p_grid2, TT, lag_max))
  } else {
    # cross correlation of the differenced process
    T0 <- object$T0
    if (hasArg(multi)) multi <- as.integer(list(...)$multi) else multi <- 1L
    ccfh <- as.numeric(ccf_sample_dp(h, x_grid1, p_grid1, x_grid2, p_grid2, T0, TT, lag_max, multi))
  }
  
  return (ccfh)
}


#' ACF of Trawl processes
#'
#' theoretical autocorrelation of a trawl process
#'
#' TODO
#'
#' CITE TODO.
#' @name ccf_trawl
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
#' @export ccf_trawl
ccf_trawl <- function(object, h, dff = 0, lag_max = 25) {
  
  trawl1 <- object$trawl[[1]]
  trawl2 <- object$trawl[[2]]
  trawl_par1 <- object$trawl_par[[1]]
  trawl_par2 <- object$trawl_par[[2]]

  if (dff < 0.5) {
    # ACF of the process itself -
    # TODO: make correct when b not equal to zero
    ccfh <- as.numeric(acf_trawl_p(h, trawl, trawl_par, lag_max))
  } else {
    # ACF of the differenced process
    # TODO
    # b <- object$b
    # acfh <- as.numeric(acf_trawl_dp(h, trawl, trawl_par, b, lag_max))
  }

  return (ccfh)
}
