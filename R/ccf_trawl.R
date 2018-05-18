#' Cross correlation function
#'
#' compute empirical or theoretical cross correlation values of either the processes
#' itself or the processes of first differences
#'
#' the object is either an output object from sim_trawl for the empirical autocorrelation
#' values, or a combination of fit_trawl and fit_levy. See the examples. When computing the sample
#' ACF, it is possible to add the argument 'multi' for reducing the 
#' estimation variance. This is only useful when considering the differenced time series,
#' since the non-differenced autocorrelations are exact.
#' 
#' The sample method works for any type of time series, regularly or irregularly spaced.
#'
#' @name ccf_trawl
#' @encoding UTF-8
#' @concept trawl
#' @param object object containing all the specifications for the process, see details
#' @param h observation frequency
#' @param method "sample" for the empirical values, something else for the theoretical values
#' @param dff number of times the series has to be differenced before computing the 
#' autocorrelation, possible values are 0 (default) and 1
#' @param lag_max number of lags to consider in each direction, default 25
#' @param \dots any other passthrough parameters
#' @author Dries Cornilly
#' @seealso \code{\link{fit_trawl}}, \code{\link{sim_trawl}}, \code{\link{acf_trawl}}
#' @references
#' Barndorff‐Nielsen, O. E., Lunde, A., Shephard, N., & Veraart, A. E. (2014). 
#' Integer‐valued Trawl Processes: A Class of Stationary Infinitely Divisible Processes. 
#' Scandinavian Journal of Statistics, 41(3), 693-724.
#' 
#' Veraart, A. E. (2018). 
#' Modelling, simulation and inference for multivariate time series of counts. 
#' arXiv preprint arXiv:1608.03154.
#'
#' @examples
#' ### multivariate negative binomial trawl process
#' trawl <- list("gamma", "exp")
#' trawl_par <- list(c(0.9, 1.8), c(0.5))
#' levy_par <- matrix(c(1.5, 0.4, 0.3,
#'                      0.7, 0.4, NA,
#'                      0.6, 0.3, NA), ncol = 3, byrow = TRUE)
#' design_matrix <- matrix(c(1, 1, 1, 0, 0, 1), nrow = 2)
# 
#' sm <- sim_trawl(list("levy_seed" = "negBin", "levy_par" = levy_par,
#'                      "trawl" = trawl, "trawl_par" = trawl_par,
#'                      "design_matrix" = design_matrix, "b" = c(0, 0),
#'                      "T0" = 0, "TT" = 25400, "observed_freq" = 1e-6),
#'                 univariate = FALSE)
#' 
#' # cross correlation
#' lag_max <- 10
#' ht <- 0.5
#' plot(seq(-lag_max * ht, lag_max * ht, ht), ccf_trawl(sm, ht, lag_max = lag_max), 
#'      pch = 16, main = "cross correlation")
#' lines(seq(-lag_max * ht, lag_max * ht, ht),
#'       ccf_trawl(list("trawl" = trawl, "trawl_par" = trawl_par, "levy_seed" = "negBin", 
#'                      "levy_par" = levy_par, "design_matrix" = design_matrix), 
#'                 h = ht, method = "theor", dff = 0, lag_max = lag_max), lwd = 2, col = "blue")
#' 
#' # cross correlation of differenced series
#' lag_max <- 10
#' ht <- 0.5
#' plot(seq(-lag_max * ht, lag_max * ht, ht), ccf_trawl(sm, ht, dff = 1, lag_max = lag_max), 
#'      pch = 16, main = "cross correlation differenced series")
#' lines(seq(-lag_max * ht, lag_max * ht, ht),
#'       ccf_trawl(list("trawl" = trawl, "trawl_par" = trawl_par,
#'                      "levy_seed" = "negBin", "levy_par" = levy_par,
#'                      "design_matrix" = design_matrix, "b" = c(0, 0)), 
#'                 h = ht, method = "theor", dff = 1, lag_max = lag_max), lwd = 2, col = "blue")
#'
#' ### multivariate Poisson / Skellam trawl process
#' trawl <- list("gamma", "exp")
#' trawl_par <- list(c(0.9, 1.8), c(0.5))
#' levy_par <- matrix(c(0.13, 0.13, 0.23, 0.11, 0.05), ncol = 1)
#' design_matrix <- matrix(c(1, 0, 0, 1, -1, 1, 1, -1, 1, 0), nrow = 2)
# 
#' sm <- sim_trawl(list("levy_seed" = "Poisson", "levy_par" = levy_par,
#'                      "trawl" = trawl, "trawl_par" = trawl_par,
#'                      "design_matrix" = design_matrix, "b" = c(0, 0), 
#'                      "T0" = 35, "TT" = 24600, "observed_freq" = 1e-6),
#'                 univariate = FALSE)
#' 
#' # cross correlation
#' lag_max <- 10
#' ht <- 0.5
#' plot(seq(-lag_max * ht, lag_max * ht, ht), ccf_trawl(sm, ht, lag_max = lag_max), 
#'      pch = 16, main = "cross correlation")
#' lines(seq(-lag_max * ht, lag_max * ht, ht),
#'       ccf_trawl(list("trawl" = trawl, "trawl_par" = trawl_par,
#'                      "levy_seed" = "Poisson", "levy_par" = levy_par,
#'                      "design_matrix" = design_matrix), 
#'                 h = ht, method = "theor", dff = 0, lag_max = lag_max), lwd = 2, col = "blue")
#' 
#' # cross correlation of differenced series
#' lag_max <- 10
#' ht <- 1
#' plot(seq(-lag_max * ht, lag_max * ht, ht), ccf_trawl(sm, ht, dff = 1, lag_max = lag_max), 
#'      pch = 16, main = "cross correlation differenced series")
#' lines(seq(-lag_max * ht, lag_max * ht, ht),
#'       ccf_trawl(list("trawl" = trawl, "trawl_par" = trawl_par,
#'                      "levy_seed" = "Poisson", "levy_par" = levy_par,
#'                      "design_matrix" = design_matrix, "b" = c(0, 0)), 
#'                 h = ht, method = "theor", dff = 1, lag_max = lag_max), lwd = 2, col = "blue")
#'
#' @export ccf_trawl
ccf_trawl <- function(object, h, method = "sample", dff = 0, lag_max = 25, ...) {
  
  if (method == "sample") {
    
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
    
  } else {
    
    # any other method means theoretical cross correlation is computed
    trawl1 <- object$trawl[[1]]
    trawl2 <- object$trawl[[2]]
    trawl_par1 <- object$trawl_par[[1]]
    trawl_par2 <- object$trawl_par[[2]]
    levy_seed <- object$levy_seed
    levy_par <- object$levy_par
    design_matrix <- object$design_matrix
    
    if (dff < 0.5) {
      # ACF of the process itself -
      # TODO: make correct when b not equal to zero
      ccfh <- as.numeric(ccf_trawl_p(h, trawl1, trawl_par1, trawl2, trawl_par2, 
                                     levy_seed, levy_par, design_matrix, lag_max))
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
      ccfh <- as.numeric(ccf_trawl_dp(h, trawl1, trawl_par1, trawl2, trawl_par2, b,
                                      levy_seed, levy_par, design_matrix, lag_max))
    }
  }
  
  return (ccfh)
}
