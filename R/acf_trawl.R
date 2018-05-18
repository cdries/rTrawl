#' Autocorrelation function
#'
#' compute empirical or theoretical autocorrelation values of either the process
#' itself or the process of first differences
#'
#' the object is either an output object from sim_trawl for the empirical autocorrelation
#' values, or an output object from fit_trawl for the theoretical ones. When computing the sample
#' ACF, it is possible to add the argument 'multi' for reducing the 
#' estimation variance. This is only useful when considering the differenced time series,
#' since the non-differenced autocorrelations are exact.
#' 
#' The sample method works for any type of time series, regularly or irregularly spaced.
#'
#' @name acf_trawl
#' @encoding UTF-8
#' @concept trawl
#' @param object object containing all the specifications for the process, see details
#' @param h observation frequency
#' @param method "sample" for the empirical values, something else for the theoretical values
#' @param dff number of times the series has to be differenced before computing the 
#' autocorrelation, possible values are 0 (default) and 1
#' @param lag_max number of lags to consider, default 25
#' @param drop_zero whether to exclude the autocorrelation at lag 0, default TRUE
#' @param \dots any other passthrough parameters
#' @author Dries Cornilly
#' @seealso \code{\link{fit_trawl}}, \code{\link{sim_trawl}}
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
#' # empirical and estimated autocorrelations
#' sim <- sim_trawl(list())
#' sim$h <- 0.5
#' sim$trawl <- "exp"
#' sim$lag_max <- 3
#' sim$method <- "acf"
#' ft <- fit_trawl(sim)
#' plot(1:10, acf_trawl(sim, sim$h, lag_max = 10))
#' lines(1:10, acf_trawl(sim, sim$h, method = "acf", lag_max = 10), col = "blue")
#' 
#' # empirical and estimated autocorrelation for a differenced process
#' sim <- sim_trawl(list("levy_seed" = "Skellam", levy_par = c(0.13, 0.11), b = 0.3))
#' sim$h <- exp(seq(log(1e-2), log(60), length.out = 51))
#' sim$trawl <- "exp"
#' sim$include_b <- TRUE
#' sim$include_cum1 <- FALSE
#' sim$method <- "vs_C"
#' ftC <- fit_trawl(sim)
#' h <- 0.5
#' plot(1:10, acf_trawl(sim, h, dff = 1, lag_max = 10))
#' lines(1:10, acf_trawl(sim, h, method = "vs_C", dff = 1, lag_max = 10), col = "blue")
#'
#' @export acf_trawl
acf_trawl <- function(object, h, method = "sample", dff = 0, lag_max = 25, drop_zero = TRUE, ...) {
  
  if (method == "sample") {
    
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
