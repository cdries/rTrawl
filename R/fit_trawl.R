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

  if (method == "vs_C") {

    # variance signature plot - cum2 and cum1 as explicit parameters
    include_cum1 <- object$include_cum1
    include_b <- object$include_b
    h <- object$h
    lfit <- fit_trawl_vs_C(h, x_grid, p_grid, T0, TT, trawl, include_cum1, include_b, ...)

  } else if (method == "vs_SY") {

    # variance signature plot - cum2 and cum1 as function of b, beta_0 and alpha_y
    include_cum1 <- object$include_cum1
    include_b <- object$include_b
    h <- object$h
    lfit <- fit_trawl_vs_SY(as.numeric(h), as.numeric(x_grid), as.numeric(p_grid), as.numeric(T0), 
                            as.numeric(TT), trawl, include_cum1, include_b, ...)

  } else if (method == "acf") {

    # if (include_b) stop("The parameter b cannot be included when estimating based on autocorrelation.")
    # if (trawl == "nonpar") {
    #   lfit  <- fit_trawl_ac_nonpar(x_grid, p_grid, T0, TT, levy_seed) # TODO
    # } else {
    #   lfit <- fit_trawl_ac_par(x_grid, p_grid, T0, TT, levy_seed, trawl, ...) # TODO
    # }

  } else {
    stop("Select valid method (vs_C, vs_SY or acf).")
  }

  return(lfit)
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


# fit_trawl_varsig_nonpar <- function(h, x_grid, p_grid, T0, TT, include_b, s0, type, ...) { # TODO: add include_xi option
# 
#   # constants
#   beta_0 <- sum(abs(diff(p_grid)) > sqrt(.Machine$double.eps)) / (TT - T0)
#   levy_alpha <- levy_inst_prob(p_grid)
#   xi1 <- beta_0 * sum(levy_alpha[, 1] * levy_alpha[, 2])
#   xi2 <- beta_0 * sum(levy_alpha[, 1]^2 * levy_alpha[, 2])
#   if (hasArg(multi)) multi <- list(...)$multi else multi <- 10
#   sig_emp <- varsig_sample(list("x_grid" = x_grid, "p_grid" = p_grid, "T0" = T0, "TT" = TT),
#                            h, multi = multi)
# 
#   # derivative of variance signature plot
#   sig2_emp <- sig_emp * h
#   n_h <- length(h)
#   d_h <- diff(h)
#   fd_sig2_emp <- rep(NA, n_h)
#   for (ii in 2:(n_h - 1)) {
#     fd_sig2_emp[ii] <- -d_h[ii] * sig2_emp[ii - 1] / (d_h[ii - 1] * (d_h[ii] + d_h[ii - 1])) +
#       (d_h[ii] - d_h[ii - 1]) * sig2_emp[ii] / (d_h[ii] * d_h[ii - 1]) +
#       d_h[ii - 1] * sig2_emp[ii + 1] / (d_h[ii] * (d_h[ii] + d_h[ii - 1]))
#   }
#   fd_sig2_emp[1] <- (sig2_emp[2] - sig2_emp[1]) / d_h[1]
#   fd_sig2_emp[n_h] <- (sig2_emp[n_h] - sig2_emp[n_h - 1]) / d_h[n_h - 1]
#   if (hasArg(nknots)) nknots <- list(...)$nknots else nknots <- max(20, round(length(h) / 3))
#   cbs_sig2_e_emp <- cobs::cobs(log(h), fd_sig2_emp - 2 * h * xi1^2, constraint = "decrease",
#                                nknots = nknots, print.mesg = FALSE, degree = 2,
#                                pointwise = matrix(c(0, -15, xi2), nrow = 1))
#   d_sig2_e <- function(x) {
#     ind <- x < .Machine$double.eps
#     val <- rep(NA, length(x))
#     val[ind] <- cobs:::.splValue(2, cbs_sig2_e_emp$knots, cbs_sig2_e_emp$coef, exp(-15), 0)
#     val[!ind] <- cobs:::.splValue(2, cbs_sig2_e_emp$knots, cbs_sig2_e_emp$coef, log(x[!ind]), 0)
#     return (val)
#   }
#   d_sig2_e_grad <- function(x) cobs:::.splValue(2, cbs_sig2_e_emp$knots, cbs_sig2_e_emp$coef, log(x), 1) / x
# 
#   # s0 and related constants
#   ind <- h < s0
#   int_d_sig2_e <- leb_NONPAR_partial(h, d_sig2_e, s0, ind)
#   int_d_sig2_e_0 <- leb_NONPAR_partial(0, d_sig2_e, s0, TRUE)
#   hs0_e <- d_sig2_e(s0)
# 
#   # boundaries
#   if (type == "exp") lb <- c(1e-6, 1e-6, -Inf) else lb <- c(1 + 1e-8, 1e-6, -Inf)
#   ub <- c(Inf, Inf, s0 - 1e-9)
#   if (include_b) {
#     lb <- c(lb, 0)
#     ub <- c(ub, 1 - 1e-7)
#   }
#   n_theta <- length(lb)
# 
#   # objective function
#   obj <- function(theta) {
# 
#     lambda <- theta[1]
#     gamma <- theta[2]
#     eta <- theta[3]
#     if (include_b) b <- theta[4] else b <- 0
# 
#     sig_theor <- varsig_trawl_nonpar_obj(lambda, gamma, eta, b, s0, beta_0, h,
#                                          int_d_sig2_e, int_d_sig2_e_0, ind, levy_alpha, type)
# 
#     # objective value
#     diffsig <- sig_theor - sig_emp
#     val <- sum(diffsig^2)
# 
#     # gradient
#     vs_grad <- varsig_trawl_nonpar_grad(lambda, gamma, eta, b, s0, beta_0, h, int_d_sig2_e,
#                                         int_d_sig2_e_0, ind, levy_alpha, type, include_b)
#     grad <- 2 * colSums(matrix(diffsig, nrow = n_h, ncol = n_theta, byrow = FALSE) * vs_grad)
# 
#     return (list("objective" = val, "gradient" = grad))
#   }
# 
#   # constraints function
#   g_eq <- function(theta) {
# 
#     lambda <- theta[1]
#     gamma <- theta[2]
#     eta <- theta[3]
#     if (include_b) b <- theta[4] else b <- 0
# 
#     # objective value
#     if (type == "exp") {
#       val <- 0.5 * ((2 - b) * hs0_e / xi2 - b) / (1 - b) - gamma * exp(-lambda * (s0 - eta))
#     } else {
#       val <- 0.5 * ((2 - b) * hs0_e / xi2 - b) / (1 - b) - gamma * (s0 - eta)^(-lambda)
#     }
# 
#     # gradient
#     grad <- rep(NA, 3 + include_b)
#     if (type == "exp") {
#       grad[1] <- gamma * (s0 - eta) * exp(-lambda * (s0 - eta))
#       grad[2] <- -exp(-lambda * (s0 - eta))
#       grad[3] <- -gamma * lambda * exp(-lambda * (s0 - eta))
#     } else {
#       grad[1] <- gamma * log(s0 - eta) * (s0 - eta)^(-lambda)
#       grad[2] <- -(s0 - eta)^(-lambda)
#       grad[3] <- -gamma * lambda * (s0 - eta)^(-lambda - 1)
#     }
#     if (include_b) grad[4] <- 0.5 * (-(1 - b) * (hs0_e / xi2 + 1) +
#                                        (2 - b) * hs0_e / xi2 - b) / (1 - b)^2
# 
#     return (list("constraints" = val, "jacobian" = grad))
#   }
# 
#   # initial estimates
#   eta0 <- 0
#   if (type == "exp") {
#     lambda0 <- 0.5
#     if (include_b) {
#       if (hasArg(b0)) b0 <- list(...)$b0 else b0 <- 0.1
#       fn_x0 <- function(x) g_eq(c(x, 1, eta0, b0))$constraints
#       sol_x0 <- nleqslv::nleqslv(x = lambda0, fn = fn_x0)
#       x0 <- c(sol_x0$x, 1, eta0, b0)
#     } else {
#       b0 <- 0
#       fn_x0 <- function(x) g_eq(c(x, 1, eta0, b0))$constraints
#       sol_x0 <- nleqslv::nleqslv(x = lambda0, fn = fn_x0)
#       x0 <- c(sol_x0$x, 1, eta0)
#     }
#   } else {
#     lambda0 <- 3
#     if (include_b) {
#       if (hasArg(b0)) b0 <- list(...)$b0 else b0 <- 0.1
#       fn_x0 <- function(x) g_eq(c(lambda0, x, eta0, b0))$constraints
#       sol_x0 <- nleqslv::nleqslv(x = 1, fn = fn_x0)
#       x0 <- c(lambda0, sol_x0$x, eta0, b0)
#     } else {
#       b0 <- 0
#       fn_x0 <- function(x) g_eq(c(x, 1, eta0, b0))$constraints
#       sol_x0 <- nleqslv::nleqslv(x = 1, fn = fn_x0)
#       x0 <- c(lambda0, sol_x0$x, eta0)
#     }
#   }
# 
#   # optimization
#   optscontrol <- list(algorithm = "NLOPT_LD_SLSQP", xtol_rel = 1e-6, xtol_abs = 1e-10,
#                       maxeval = 10000, print_level = 0, check_derivatives = FALSE)
#   sol <- nloptr::nloptr(x0 = x0, eval_f = obj, lb = lb, ub = ub, eval_g_eq = g_eq, opts = optscontrol)
#   lambda <- sol$solution[1]
#   gamma <- sol$solution[2]
#   eta <- sol$solution[3]
#   if (include_b) b <- sol$solution[4] else b <- 0
# 
#   # trawl function
#   d_trawl <- function(x) {
#     d <- rep(0, length(x))
#     ind <- which((x > .Machine$double.eps) * (x < s0 + .Machine$double.eps) > 0.5)
#     d[which((x > -.Machine$double.eps) * (x < .Machine$double.eps) > 0.5)] <- 1
#     d[ind] <- b + 0.5 * ((2 - b) * d_sig2_e(x[ind]) / xi2 - b)
#     if (type == "exp") {
#       d[x > s0] <- b + (1 - b) * gamma * exp(-lambda * (x[x > s0] - eta))
#     } else {
#       d[x > s0] <- b + (1 - b) * gamma * (x[x > s0] - eta)^(-lambda)
#     }
# 
#     return (d)
#   }
# 
#   # return objects
#   trawl_par <- list("d_trawl" = d_trawl, "lambda" = lambda, "gamma" = gamma,
#                     "eta" = eta, "s0" = s0, "type" = type)
#   levy_par <- levy_alpha2nu(levy_alpha, b, beta_0)
# 
#   return (list("trawl" = "nonpar", "trawl_par" = trawl_par, "b" = b,
#                "levy_seed" = "nonpar", "levy_par" = levy_par))
# }
# 
# 
# fit_trawl_racf <- function(h, x_grid, p_grid, T0, TT, trawl, include_b, include_xi, ...) {
# 
#   # contants
#   n_trawl <- number_parameters_trawl(trawl)
#   n_h <- length(h)
#   if (hasArg(multi)) multi <- list(...)$multi else multi <- 10
#   ar1_emp <- rep(0, n_h)
#   for (ii in 1:n_h) ar1_emp[ii] <- acf_sample(list("x_grid" = x_grid, "p_grid" = p_grid, "T0" = T0, "TT" = TT),
#                                               h[ii], dff = 1, lag_max = 1, drop_zero = TRUE, multi = multi)
# 
#   # bounds
#   bounds <- trawl_bounds(trawl)
#   lb <- bounds$lb
#   ub <- bounds$ub
#   if (include_b) {
#     lb <- c(lb, 0)
#     ub <- c(ub, 1)
#   }
#   if (include_xi) {
#     lb <- c(lb, 0)
#     ub <- c(ub, Inf)
#   }
#   n_theta <- length(lb)
# 
#   # initial estimate
#   if (hasArg(x0)) {
#     x0 <- list(...)$x0
#   } else {
#     x0 <- trawl_x0(trawl)
#     if (include_b) x0 <- c(x0, 0.5)
#     if (include_xi) x0 <- c(x0, 0)
#   }
# 
#   # objective function
#   obj <- function(theta) {
# 
#     trawl_par <- theta[1:n_trawl]
#     if (include_b) b <- theta[n_trawl + 1] else b <- 0
#     if (include_xi) xi <- theta[n_theta] else xi <- 0
# 
#     # objective value
#     temp <- ar1_trawl_obj_grad(trawl, trawl_par, b, h, include_b, include_xi, xi)
#     ar1_theor <- temp$ar1
#     diffar1 <- ar1_theor - ar1_emp
#     val <- sum(diffar1^2)
# 
#     # gradient
#     ar1_grad <- temp$ar1_grad
#     grad <- 2 * colSums(matrix(diffar1, nrow = n_h, ncol = n_theta, byrow = FALSE) * ar1_grad)
# 
#     return (list("objective" = val, "gradient" = grad))
#   }
# 
#   # optimization
#   optscontrol <- list(algorithm = "NLOPT_LD_MMA", xtol_rel = 1e-05, maxeval = 10000, print_level = 0,
#                       check_derivatives = FALSE)
#   sol <- nloptr::nloptr(x0 = x0, eval_f = obj, lb = lb, ub = ub, opts = optscontrol)
# 
#   # return object
#   trawl_par <- sol$solution[1:n_trawl]
#   if (include_b) b <- sol$solution[n_trawl + 1] else b <- 0
# 
#   return (list("trawl" = trawl, "trawl_par" = trawl_par, "b" = b))
# }
# 
# 
# fit_trawl_vs <- function(h, x_grid, p_grid, T0, TT, trawl, ...) {
# 
#   # contants
#   n_trawl <- number_parameters_trawl(trawl)
#   n_h <- length(h)
#   if (hasArg(multi)) multi <- list(...)$multi else multi <- 5
#   sig_emp <- vs_sample(x_grid, p_grid, T0, TT, h, multi = multi)
# 
#   # bounds
#   bounds <- trawl_bounds(trawl)
#   lb <- c(bounds$lb, 0, 0)
#   ub <- c(bounds$ub, Inf, Inf)
#   n_theta <- length(lb)
# 
#   # initial estimate
#   if (hasArg(x0)) {
#     x0 <- list(...)$x0
#   } else {
#     x0 <- c(trawl_x0(trawl), 0.01, 0.01)
#   }
# 
#   # objective function
#   obj <- function(theta) {
# 
#     trawl_par <- theta[1:n_trawl]
#     omega <- theta[n_trawl + 1]
#     xi <- theta[n_trawl + 2]
# 
#     # objective value
#     temp <- vs_obj_grad(trawl, trawl_par, omega, xi, h)
#     diffsig <- temp$sig_theor - sig_emp
#     val <- sum(diffsig^2)
# 
#     # gradient
#     vs_grad <- temp$vs_grad
#     grad <- 2 * colSums(matrix(diffsig, nrow = n_h, ncol = n_theta, byrow = FALSE) * vs_grad)
# 
#     return (list("objective" = val, "gradient" = grad))
#   }
# 
#   # optimization
#   optscontrol <- list(algorithm = "NLOPT_LD_MMA", xtol_rel = 1e-05, maxeval = 10000,
#                       print_level = 0, check_derivatives = FALSE)
#   sol <- nloptr::nloptr(x0 = x0, eval_f = obj, lb = lb, ub = ub, opts = optscontrol)
# 
#   # return object
#   trawl_par <- sol$solution[1:n_trawl]
#   omega <- sol$solution[n_trawl + 1]
#   xi <- sol$solution[n_trawl + 2]
# 
#   return (list("trawl" = trawl, "trawl_par" = trawl_par, "omega" = omega, "xi" = xi, "sol" = sol))
# }
# 
