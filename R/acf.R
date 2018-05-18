#' Autocorrelation function
#'
#' theoretical autocorrelation of a trawl process
#'
#' TODO
#'
#' CITE TODO.
#' @name acf_trawl
#' @concept trawl
#' @param object bla
#' @param h bla
#' @param method bla
#' @param dff bla
#' @param lag_max bla
#' @param drop_zero bla
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
#' @export acf_trawl
acf_trawl <- function(object, h, method = "sample", dff = 0, lag_max = 25, drop_zero = TRUE) {
  
  if ("method" == sample) {
    
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
    
  } else {
    
    # any other method means theoretical autocorrelation
    trawl <- object$trawl
    trawl_par <- object$trawl_par
    
    if (dff < 0.5) {
      # ACF of the process itself - 
      # TODO: make correct when b not equal to zero
      acfh <- as.numeric(acf_trawl_p(h, trawl, trawl_par, lag_max))
    } else {
      # ACF of the differenced process - compute b when estimation method was "vs_C"
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
  }
  
  if (!drop_zero) acfh <- c(1, acfh)
  
  return (acfh)
}
