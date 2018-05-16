#' Sample ACF of irregularly observed time series
#'
#' sample autocorrelation of a trawl process
#'
#' TODO
#'
#' CITE TODO.
#' @name acf_sample
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
#' @export acf_sample
#' @useDynLib rTrawl
acf_sample <- function(object, h, dff = 0, lag_max = 25, drop_zero = TRUE, ...) {
  
  x_grid <- object$x_grid
  p_grid <- object$p_grid
  TT <- object$TT
  
  if (dff < 0.5) {
    # ACF of the process itself
    acfh <- as.numeric(acf_sample_p(h, x_grid, p_grid, TT, lag_max))
  } else {
    # ACF of the differenced process
    T0 <- object$T0
    if (hasArg(multi)) multi <- as.integer(list(...)$multi) else multi <- 1L
    acfh <- as.numeric(acf_sample_dp(h, x_grid, p_grid, T0, TT, lag_max, multi))
  }
  
  if (!drop_zero) acfh <- c(1, acfh)
  
  return (acfh)
}


#' ACF of Trawl processes
#'
#' theoretical autocorrelation of a trawl process
#'
#' TODO
#'
#' CITE TODO.
#' @name acf_trawl
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
#' @export acf_trawl
acf_trawl <- function(object, h, dff = 0, lag_max = 25, drop_zero = TRUE) {
  
  trawl <- object$trawl
  trawl_par <- object$trawl_par
  
  if (dff < 0.5) {
    # ACF of the process itself - 
    # TODO: make correct when b not equal to zero
    acfh <- as.numeric(acf_trawl_p(h, trawl, trawl_par, lag_max))
  } else {
    # ACF of the differenced process
    b <- object$b
    if (is.null(b)) {
      if (is.null(object$xi)) {
        b <- 0
      } else {
        b <- object$xi / object$omega
        b <- b / (1 + b)
      }
    }
    acfh <- as.numeric(acf_trawl_dp(h, trawl, trawl_par, b, lag_max))
  }
  
  if (!drop_zero) acfh <- c(1, acfh)
  
  return (acfh)
}

