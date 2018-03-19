#' Estimation of Trawl processes
#'
#' estimates a path of a Trawl process
#'
#' TODO
#'
#' CITE TODO.
#' @name fit_trawl
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
#' @export fit_trawl
fit_trawl <- function(object, ...) {
  
  # extract settings
  method <- object$method
  trawl <- object$trawl
  x_grid <- object$x_grid
  p_grid <- object$p_grid
  T0 <- object$T0
  TT <- object$TT
  h <- object$h
  
  if (method == "vs_C") {
    
    # variance signature plot - cum2_levy, cum2_trawl and cum1_levy as explicit parameters
    include_cum1 <- object$include_cum1
    include_b <- object$include_b # if false, only cum2_trawl is a parameter
    lfit <- fit_trawl_vs_C(as.numeric(h), as.numeric(x_grid), as.numeric(p_grid), as.numeric(T0), 
                           as.numeric(TT), trawl, include_cum1, include_b, ...)
    
  } else if (method == "vs_SY") {
    
    # variance signature plot - cum2_Levy and cum1 as function of b, beta_0 and alpha_y
    include_cum1 <- object$include_cum1
    include_b <- object$include_b
    lfit <- fit_trawl_vs_SY(as.numeric(h), as.numeric(x_grid), as.numeric(p_grid), as.numeric(T0), 
                            as.numeric(TT), trawl, include_cum1, include_b, ...)
    
  } else if (method == "acf") {
    
    # first fit the trawl
    # TODO - add check for multivariate
    # TODO - allow multiple h values for each trawl when fitting multivariate
    # TODO - multivariate is just looping over the univariate cases
    lag_max <- object$lag_max
    lfit <- fit_trawl_acf(as.numeric(h), as.integer(lag_max), as.numeric(x_grid), as.numeric(p_grid), 
                          as.numeric(T0), as.numeric(TT), trawl, ...)
    
  } else {
    stop("Select valid method (vs_C, vs_SY or acf).")
  }
  
  return(lfit)
}


#' @import nloptr
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


#' @import nloptr
fit_trawl_vs_SY <- function(h, x_grid, p_grid, T0, TT, trawl, include_cum1, include_b, ...) {
  
  # contants
  n_trawl <- number_parameters_trawl(trawl)
  n_h <- length(h)
  if (hasArg(multi)) multi <- as.integer(list(...)$multi) else multi <- 1L
  vs_emp <- vs_sample(h, x_grid, p_grid, T0, TT, multi)
  tmp <- levy_alpha_beta(p_grid, T0, TT)
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
               "levy_seed" = "nonpar", "levy_par" = levy_par))
}


#' @import nloptr
fit_trawl_acf <- function(h, lag_max, x_grid, p_grid, T0, TT, trawl, ...) {
  
  # contants
  n_trawl <- number_parameters_trawl(trawl)
  if (hasArg(multi)) multi <- as.integer(list(...)$multi) else multi <- 1L
  acf_emp <- acf_sample_p(h, x_grid, p_grid, T0, TT, lag_max, multi)
  
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

  return (list("trawl" = trawl, "trawl_par" = trawl_par))
}


