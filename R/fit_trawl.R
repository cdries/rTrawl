#' Estimation of Trawl processes
#'
#' fits chosen trawl to a given univariate trawl process, multivariate process should be passed
#' individually sequantially.
#'
#' the passed object contains all the process specifications in list format. The following
#' arguments should ideally be present in the object: vector with time points of the observations
#' (x_grid), vector with process values (p_grid), trawl (trawl), corresponding trawl parameters (trawl_par),
#' initial observation time (T0), end of observation period (TT), vector with observation
#' frequencies (h) when fitting with "vs_C" or "vs_SY", or single observation frequency when fitting
#' with "acf". In addition, the estimation method should be passed, "vs_C" (default), "vs_SY" or "acf",
#' we refer to the literature for further details on these.
#' 
#' When estamating based on the variance signature plot ("vs_C" or "vs_SY"), it is also
#' possible to specify whether the first cumulant (include_cum1") should be included (default FALSE),
#' and whether or not a pure levy component is present (include b, default TRUE).
#' 
#' When estimating based on the autocorrelation function, there should not be a pure levy
#' component present, and it is possible to set the maximum number of lag used by
#' "lag_max", default equals 5.
#' 
#' All estimation methods accept the passthrough arguument 'multi' for reducing the 
#' estimation variance.
#' 
#' See examples for the three methods of estimation.
#'
#' @name fit_trawl
#' @concept trawl
#' @encoding UTF-8
#' @param object object containing all the specifications for the process, see details
#' @param \dots any other passthrough parameters
#' @author Dries Cornilly
#' @seealso \code{\link{sim_trawl}}, \code{\link{vs_trawl}}, \code{\link{acf_trawl}}
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
#' # estimation based on the variance signature plot
#' sim <- sim_trawl(list("levy_seed" = "Skellam", levy_par = c(0.13, 0.11), b = 0.3))
#' sim$h <- exp(seq(log(1e-2), log(60), length.out = 51))
#' sim$trawl <- "exp"
#' sim$include_b <- TRUE
#' sim$include_cum1 <- FALSE
#' 
#' sim$method <- "vs_C"
#' ftC <- fit_trawl(sim)
#' 
#' sim$method <- "vs_SY"
#' ftSY <- fit_trawl(sim)
#' 
#' plot(log(sim$h), vs_trawl(sim, sim$h))
#' lines(log(sim$h), vs_trawl(ftC, sim$h, method = "vs_C"), col = "blue")
#' lines(log(sim$h), vs_trawl(ftSY, sim$h, method = "vs_SY"), col = "red")
#' legend("topright", c("vs_C", "vs_SY"), lwd = c(2, 2), col = c("blue", "red"))
#' 
#' # estimation based on the autocorrelation
#' sim <- sim_trawl(list())
#' sim$h <- 0.5
#' sim$trawl <- "exp"
#' sim$lag_max <- 3
#' sim$method <- "acf"
#' ft <- fit_trawl(sim)
#' plot(1:10, acf_trawl(sim, sim$h, lag_max = 10))
#' lines(1:10, acf_trawl(ft, sim$h, method = "acf", lag_max = 10), col = "blue")
#'
#' @export fit_trawl
fit_trawl <- function(object, ...) {
  
  # extract settings
  nn <- names(object)
  x_grid <- object$x_grid
  p_grid <- object$p_grid
  if ("method" %in% nn) method <- object$method else method <- "vs_C"
  if ("trawl" %in% nn) trawl <- object$trawl else trawl <- "exp"
  if ("T0" %in% nn) T0 <- object$T0 else T0 <- 0
  if ("TT" %in% nn) TT <- object$TT else TT <- 3600
  if ("h" %in% nn) h <- object$h else h <- exp(seq(log(1e-2), log(60), length.out = 51))

  if (method == "vs_C") {
    
    # variance signature plot - cum2_levy, cum2_trawl and cum1_levy as explicit parameters
    if ("include_cum1" %in% nn) include_cum1 <- object$include_cum1 else include_cum1 <- FALSE
    if ("include_b" %in% nn) include_b <- object$include_b else include_b <- TRUE # if false, only cum2_trawl is a parameter
    lfit <- fit_trawl_vs_C(as.numeric(h), as.numeric(x_grid), as.numeric(p_grid), as.numeric(T0), 
                           as.numeric(TT), trawl, include_cum1, include_b, ...)
    
  } else if (method == "vs_SY") {
    
    # variance signature plot - cum2_Levy and cum1 as function of b, beta_0 and alpha_y
    if ("include_cum1" %in% nn) include_cum1 <- object$include_cum1 else include_cum1 <- FALSE
    if ("include_b" %in% nn) include_b <- object$include_b else include_b <- TRUE
    lfit <- fit_trawl_vs_SY(as.numeric(h), as.numeric(x_grid), as.numeric(p_grid), as.numeric(T0), 
                            as.numeric(TT), trawl, include_cum1, include_b, ...)
    
  } else if (method == "acf") {
    
    # first fit the trawl
    # TODO - add check for multivariate
    # TODO - allow multiple h values for each trawl when fitting multivariate
    # TODO - multivariate is just looping over the univariate cases
    if ("lag_max" %in% nn) lag_max <- object$lag_max else lag_max <- 5
    lfit <- fit_trawl_acf(as.numeric(h), as.integer(lag_max), as.numeric(x_grid), as.numeric(p_grid), 
                          as.numeric(T0), as.numeric(TT), trawl, ...)
    
  } else {
    stop("Select valid method (vs_C, vs_SY or acf).")
  }
  
  return(lfit)
}


#' @importFrom nloptr nloptr
fit_trawl_vs_C <- function(h, x_grid, p_grid, T0, TT, trawl, include_cum1, include_b, ...) {
  
  # contants
  n_trawl <- number_parameters_trawl(trawl)
  n_h <- length(h)
  if (hasArg(multi)) multi <- as.integer(list(...)$multi) else multi <- 1L
  vs_emp <- vs_sample(h, x_grid, p_grid, T0, TT, multi)
  
  # bounds
  bounds <- trawl_bounds(trawl)
  lb <- c(bounds$lb, 0)
  ub <- c(bounds$ub, Inf)
  if (include_b) {
    lb <- c(lb, 0)
    ub <- c(ub, Inf)
    if (include_cum1) {
      lb <- c(lb, 0)
      ub <- c(ub, Inf)
    }
  }
  n_theta <- length(lb)
  
  # initial estimate
  if (hasArg(x0)) {
    x0 <- list(...)$x0
  } else {
    x0 <- c(as.numeric(trawl_x0(trawl)), 0.01)
    if (include_b) x0 <- c(x0, 0.01)
    if (include_cum1) x0 <- c(x0, 1e-6)
  }
  
  # objective function
  obj <- function(theta) {
    
    trawl_par <- theta[1:n_trawl]
    omega <- theta[n_trawl + 1]
    if (include_b) xi <- theta[n_trawl + 2] else xi <- 0
    if (include_cum1) eta <- theta[n_trawl + 3] else eta <- 0
    
    # objective value
    tmp <- vs_C(h, trawl, trawl_par, omega, xi, include_b, eta, include_cum1)
    vs_diff <- tmp$vs_theor - vs_emp
    val <- sum(vs_diff^2)
    
    # gradient
    vs_grad <- tmp$vs_grad
    grad <- 2 * colSums(matrix(vs_diff, nrow = n_h, ncol = n_theta, byrow = FALSE) * vs_grad)
    
    return (list("objective" = val, "gradient" = grad))
  }
  
  # optimization
  optscontrol <- list(algorithm = "NLOPT_LD_MMA", xtol_rel = 1e-05, maxeval = 50000, 
                      print_level = 0, check_derivatives = FALSE)
  sol <- nloptr::nloptr(x0 = x0, eval_f = obj, lb = lb, ub = ub, opts = optscontrol)
  
  # return object
  trawl_par <- sol$solution[1:n_trawl]
  omega <- sol$solution[n_trawl + 1]
  if (include_b) xi <- sol$solution[n_trawl + 2] else xi <- 0
  if (include_cum1) eta <- sol$solution[n_trawl + 3] else eta <- 0
  
  return (list("trawl" = trawl, "trawl_par" = trawl_par, "omega" = omega,
               "xi" = xi, "eta" = eta))
}


#' @importFrom nloptr nloptr
fit_trawl_vs_SY <- function(h, x_grid, p_grid, T0, TT, trawl, include_cum1, include_b, ...) {
  
  # contants
  n_trawl <- number_parameters_trawl(trawl)
  n_h <- length(h)
  if (hasArg(multi)) multi <- as.integer(list(...)$multi) else multi <- 1L
  vs_emp <- vs_sample(h, x_grid, p_grid, T0, TT, multi)
  tmp <- levy_alpha_beta(x_grid, p_grid, T0, TT)
  beta_0 <- tmp$beta_0
  levy_alpha <- tmp$levy_alpha
  
  # bounds
  bounds <- trawl_bounds(trawl)
  lb <- bounds$lb
  ub <- bounds$ub
  if (include_b) {
    lb <- c(lb, 0)
    ub <- c(ub, 1)
  }
  n_theta <- length(lb)
  
  # initial estimate
  if (hasArg(x0)) {
    x0 <- list(...)$x0
  } else {
    x0 <- as.numeric(trawl_x0(trawl))
    if (include_b) x0 <- c(x0, 0.5)
  }
  
  # objective function
  obj <- function(theta) {
    
    trawl_par <- theta[1:n_trawl]
    if (include_b) b <- theta[n_trawl + 1] else b <- 0
    
    # objective value
    tmp <- vs_SY(h, trawl, trawl_par, beta_0, levy_alpha, include_cum1, b, include_b)
    vs_diff <- tmp$vs_theor - vs_emp
    val <- sum(vs_diff^2)
    
    # gradient
    vs_grad <- tmp$vs_grad
    grad <- 2 * colSums(matrix(vs_diff, nrow = n_h, ncol = n_theta, byrow = FALSE) * vs_grad)
    
    return (list("objective" = val, "gradient" = grad))
  }
  
  # optimization
  optscontrol <- list(algorithm = "NLOPT_LD_MMA", xtol_rel = 1e-05, maxeval = 10000, 
                      print_level = 0, check_derivatives = FALSE)
  sol <- nloptr::nloptr(x0 = x0, eval_f = obj, lb = lb, ub = ub, opts = optscontrol)
  
  # return object
  trawl_par <- sol$solution[1:n_trawl]
  if (include_b) b <- sol$solution[n_trawl + 1] else b <- 0
  levy_par <- levy_alpha2nu(levy_alpha, b, beta_0)
  
  return (list("trawl" = trawl, "trawl_par" = trawl_par, "b" = b,
               "levy_seed" = "nonpar", "levy_par" = levy_par, 
               "beta_0" = beta_0, "levy_alpha" = levy_alpha, "include_cum1" = include_cum1))
}


#' @importFrom nloptr nloptr
fit_trawl_acf <- function(h, lag_max, x_grid, p_grid, T0, TT, trawl, ...) {
  
  # contants
  n_trawl <- number_parameters_trawl(trawl)
  acf_emp <- acf_sample_p(h, x_grid, p_grid, TT, lag_max)
  levy_ab <- levy_alpha_beta(x_grid, p_grid, T0, TT)
  
  # bounds
  bounds <- trawl_bounds(trawl)
  lb <- bounds$lb
  ub <- bounds$ub
  n_theta <- length(lb)
  
  # initial estimate
  if (hasArg(x0)) x0 <- list(...)$x0 else x0 <- as.numeric(trawl_x0(trawl))
  
  # objective function
  obj <- function(theta) {
    
    # objective value
    tmp <- acf_BN_V(h, trawl, theta, lag_max)
    acf_diff <- tmp$acf_theor - acf_emp
    val <- sum(acf_diff^2)

    # gradient
    acf_grad <- tmp$acf_grad
    grad <- 2 * colSums(matrix(acf_diff, nrow = lag_max, ncol = n_theta, byrow = FALSE) * acf_grad)

    return (list("objective" = val, "gradient" = grad))
  }
  
  # optimization
  optscontrol <- list(algorithm = "NLOPT_LD_MMA", xtol_rel = 1e-05, maxeval = 10000, 
                      print_level = 0, check_derivatives = FALSE)
  sol <- nloptr::nloptr(x0 = x0, eval_f = obj, lb = lb, ub = ub, opts = optscontrol)
  
  # return object
  trawl_par <- sol$solution

  return (list("trawl" = trawl, "trawl_par" = trawl_par,
               "beta_0" = levy_ab$beta_0, "levy_alpha" = levy_ab$levy_alpha))
}
