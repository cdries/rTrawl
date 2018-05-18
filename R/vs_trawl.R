#' Variance signature values
#'
#' compute empirical and theoretical variance signature values
#'
#' the object is either an output object from sim_trawl for the empirical variance signature
#' values, or an output object from fit_trawl for the theoretical ones. When computing the sample
#' variance signature plot, it is possible to add the argument 'multi' for reducing the 
#' estimation variance. See Shephard and Yang (2017) for a definition of the variance 
#' signature values.
#'
#' @name vs_trawl
#' @concept trawl
#' @param object object containing all the specifications for the process, see details
#' @param h vector containing the observation frequencies at 
#' which to compute the variance signature values
#' @param method "sample" for the empirical values, "vs_C" or "vs_SY" for theoretical values when
#' estimated by these respective mtehods with fit_trawl
#' @param \dots any other passthrough parameters
#' @author Dries Cornilly
#' @seealso \code{\link{fit_trawl}}, \code{\link{sim_trawl}}
#' @references
#' Shephard, N., & Yang, J. J. (2017). 
#' Continuous time analysis of fleeting discrete price moves. 
#' Journal of the American Statistical Association, 112(519), 1090-1106.
#'
#' @examples
#' # simulate a trawl process
#' sim <- sim_trawl(list("levy_seed" = "Skellam", levy_par = c(0.13, 0.11), b = 0.3))
#' h <- exp(seq(log(1e-2), log(60), length.out = 51))
#'
#' # empirical variance signature values
#' vsS <- vs_trawl(sim, h)
#'
#' # theoretical values when estimated with "vs_C"
#' sim$method <- "vs_C"     # add fitting method
#' sim$h <- h               # add grid on which to fit
#' ft <- fit_trawl(sim)     # fit the trawl process
#' vsC <- vs_trawl(ft, h, method = "vs_C")  # compute fitted variance signature values
#'
#' # theoretical values when estimated with "vs_SY"
#' sim$method <- "vs_SY"    # add fitting method
#' sim$h <- h               # add grid on which to fit
#' ft <- fit_trawl(sim)     # fit the trawl process
#' vsSY <- vs_trawl(ft, h, method = "vs_SY")  # compute fitted variance signature values
#' 
#' # plot resulting fits
#' plot(log(h), vsS)
#' lines(log(h), vsC, col = "blue")
#' lines(log(h), vsSY, col = "red")
#' legend("topright", c("vs_C", "vs_SY"), lwd = c(2, 2), col = c("blue", "red"))
#'
#' @export vs_trawl
vs_trawl <- function(object, h, method = "sample", ...) {
  
  if (method == "sample") {

    if (hasArg(multi)) multi <- as.integer(list(...)$multi) else multi <- 1L
    vs <- vs_sample(as.numeric(h), as.numeric(object$x_grid), as.numeric(object$p_grid),
                    as.numeric(object$T0), as.numeric(object$TT), as.integer(multi))

  } else if (method == "vs_C") {

    vs <- vs_C(as.numeric(h), object$trawl, as.numeric(object$trawl_par), as.numeric(object$omega),
               as.numeric(object$xi), FALSE, as.numeric(object$eta), FALSE)$vs_theor

  } else if (method == "vs_SY") {

    vs <- vs_SY(as.numeric(h), object$trawl, as.numeric(object$trawl_par),
                as.numeric(object$beta_0), object$levy_alpha, object$include_cum1, 
                as.numeric(object$b), FALSE)$vs_theor

  } else {
    stop("provide a valid method.")
  }

  return(vs)
}
