#' Estimation of Levy basis
#'
#' estimates the parameters of a univariate or multivariate levy basis using
#' the method of moments
#'
#' the passed object contains all the process specifications in list format. The following
#' arguments should be present in the object: levy basis wanted to fit (levy_seed). Both the
#' first order moments (k1_sample) and the second order cumulants (k2_sample) may be provided.
#' If these are not present, then the sample moments are estimated and used. In order
#' to estimate the sample moments, the fitted trawls (trawl) and corresponding trawl
#' parameters (trawl_par) should be given, along with the vector with time stamps (x_grid),
#' vector with process values (p_grid), initial observation time (T0), end of observation period (TT).
#' 
#' When fitting a multivariate levy basis, the design matrix (design_matrix) should also
#' be provided, and a constraints matrix may be provided but is not required.
#' 
#' See the examples for complete use cases.
#'
#' @name fit_levy
#' @concept trawl
#' @encoding UTF-8
#' @param object object containing all the specifications for the process, see details
#' @param \dots any other passthrough parameters
#' @author Dries Cornilly
#' @seealso \code{\link{fit_trawl}}, \code{\link{sim_trawl}}
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
#' ### univariate
#' sim <- sim_trawl(list("levy_seed" = "Skellam", "levy_par" = c(0.131, 0.130),
#'                       "trawl" = "gig", "trawl_par" = c(0.01, 0.45, -0.6),
#'                       "T0" = 72.03, "TT" = 75600, "observed_freq" = 1e-6, 
#'                       "b" = 0))
#'                       
#' sim$h <- 0.5
#' sim$trawl <- "gig"
#' sim$lag_max <- 3
#' sim$method <- "acf"
#' ft <- fit_trawl(sim) 
#' 
#' lv_fit <- fit_levy(list("levy_seed" = "Skellam", "trawl" = "gig",
#'                         "trawl_par" = ft$trawl_par, "T0" = 72.03, "TT" = 75600,
#'                         "x_grid" = sim$x_grid, "p_grid" = sim$p_grid))                    
#' 
#' ### multivariate
#' ## negative binomial
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
#' h <- 1       
#' lag_max <- 2       
#' ft1 <- fit_trawl(list("T0" = 0, "TT" = 25400, "method" = "acf", "h" = h, "trawl" = trawl[[1]],
#'                       "x_grid" = sm$x_grid[[1]], "p_grid" = sm$p_grid[[1]],
#'                       "lag_max" = lag_max, "b" = 0), multi = 1)
#' ft2 <- fit_trawl(list("T0" = 0, "TT" = 25400, "method" = "acf", "h" = h, "trawl" = trawl[[2]],
#'                       "x_grid" = sm$x_grid[[2]], "p_grid" = sm$p_grid[[2]],
#'                       "lag_max" = lag_max, "b" = 0), multi = 1)
#'                       
#' lv_fit <- fit_levy(list("levy_seed" = "negBin", "trawl" = trawl,
#'                         "trawl_par" = list(ft1$trawl_par, ft2$trawl_par), "T0" = 0, "TT" = 25400,
#'                         "x_grid" = sm$x_grid, "p_grid" = sm$p_grid, 
#'                         "design_matrix" = design_matrix))
#' 
#' ## Poisson / Skellam
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
#' h <- 1
#' lag_max <- 2              
#' ft1 <- fit_trawl(list("T0" = 35, "TT" = 24600, "method" = "acf", "h" = h, "trawl" = trawl[[1]],
#'                       "x_grid" = sm$x_grid[[1]], "p_grid" = sm$p_grid[[1]],
#'                       "lag_max" = lag_max, "b" = 0), multi = 1)
#' ft2 <- fit_trawl(list("T0" = 35, "TT" = 24600, "method" = "acf", "h" = h, "trawl" = trawl[[2]],
#'                       "x_grid" = sm$x_grid[[2]], "p_grid" = sm$p_grid[[2]],
#'                       "lag_max" = lag_max, "b" = 0), multi = 1)
#' 
#' lv_fit <- fit_levy(list("levy_seed" = "Poisson", "trawl" = trawl, 
#'                         "trawl_par" = list(ft1$trawl_par, ft2$trawl_par), "T0" = 35, "TT" = 24600,
#'                         "x_grid" = sm$x_grid, "p_grid" = sm$p_grid, 
#'                         "design_matrix" = design_matrix))
#'
#' @export fit_levy
fit_levy <- function(object, ...) {
  
  # extract settings
  levy_seed <- object$levy_seed
  
  # TODO cumulants - compute if required - put in function!
  if ("k1_sample" %in% names(object)) {
    k1_sample <- object$k1_sample
    k2_sample <- object$k2_sample
  } else {
    trawl <- object$trawl
    trawl_par <- object$trawl_par
    
    if (length(trawl) == 1) {
      lebA <- leb_AtA(0.0, trawl, trawl_par)
      k1_sample <- cum_sample(1L, as.numeric(object$x_grid), as.numeric(object$p_grid), as.numeric(object$TT)) / lebA
      k2_sample <- cum_sample(2L, as.numeric(object$x_grid), as.numeric(object$p_grid), as.numeric(object$TT)) / lebA
    } else {
      k_temp <- levy_cum_mv2fit(object$T0, object$TT, object$x_grid, object$p_grid, trawl, trawl_par, length(trawl))
      k1_sample <- k_temp$k1_sample
      k2_sample <- k_temp$k2_sample
    }
  }
  
  # fit using 1st (and 2nd) cumulant
  if (length(k1_sample) == 1) {
    lfit <- levy_cum_fit(levy_seed, k1_sample, k2_sample)
  } else {
    design_matrix <- object$design_matrix
    if ("constraints" %in% names(object)) constraints <- object$constraints else constraints <- NULL
    lfit <- levy_cum_fit_mv(levy_seed, k1_sample, k2_sample, design_matrix, constraints, ...)
  }
  
  return(lfit)
}


#' @importFrom nloptr nloptr
levy_cum_fit_mv <- function(levy_seed, k1_sample, k2_sample, design_matrix, constraints, ...) {
  
  p <- nrow(design_matrix)
  p2 <- p * (p + 1) / 2
  
  if (levy_seed == "Poisson" || levy_seed == "Skellam") {
    
    if (is.null(constraints)) constraints <- diag(ncol(design_matrix))
    k <- ncol(constraints)
    
    obj <- function(theta) {
      
      theta <- constraints %*% theta
      
      k1_theor <- design_matrix %*% theta
      k2_theor <- (design_matrix %*% diag(as.numeric(theta)) %*% t(design_matrix))[lower.tri(diag(p), diag = TRUE)]
      
      k1_diff <- k1_theor - k1_sample
      k2_diff <- k2_theor - k2_sample
      
      val <- sum(k1_diff^2) + sum(k2_diff^2)
      
      gradmat <- matrix(NA, nrow = p2, ncol = k)
      kk <- 1
      for (ii in 1:p) {
        for (jj in ii:p) {
          gradmat[kk,] <- design_matrix[ii,] * design_matrix[jj,]
          kk <- kk + 1
        }
      }
      
      grad <- 2 * colSums(matrix(k1_diff, nrow = p, ncol = k, byrow = FALSE) * design_matrix) +
        2 * colSums(matrix(k2_diff, nrow = p2, ncol = k, byrow = FALSE) * gradmat)
      
      return (list("objective" = val, "gradient" = grad))
    }
    
    if (hasArg(x0)) x0 <- list(...)$x0 else x0 <- rep(0.1, k)
    optscontrol <- list(algorithm = "NLOPT_LD_MMA", xtol_rel = 1e-05, maxeval = 10000, 
                        print_level = 0, check_derivatives = FALSE)
    sol <- nloptr::nloptr(x0 = x0, eval_f = obj, lb = rep(0, k), ub = rep(Inf, k), opts = optscontrol)
    levy_par <- matrix(constraints %*% sol$solution, ncol = 1)
    
  } else if (levy_seed == "negBin") {
    
    k <- nrow(design_matrix) + ncol(design_matrix)
    
    obj <- function(theta) {
      
      theta_par <- theta2theta_par(theta, design_matrix)
      kappa_marg <- design_matrix %*% theta_par[, 1]
      alpha_marg <- theta[1:p]
      
      k1_theor <- kappa_marg * alpha_marg
      k2_theor <- levy_varcovar("negBin", theta_par, design_matrix)[lower.tri(diag(p), diag = TRUE)]
      
      k1_diff <- k1_theor - k1_sample
      k2_diff <- k2_theor - k2_sample
      
      val <- sum(k1_diff^2) + sum(k2_diff^2)
      
      k1_grad <- cbind(diag(c(kappa_marg)), design_matrix * 
                         matrix(alpha_marg, nrow = nrow(design_matrix), ncol = ncol(design_matrix), byrow = FALSE))
      k2_grad <- matrix(0, nrow = p2, ncol = k)
      kk <- 1
      for (ii in 1:p) {
        for (jj in ii:p) {
          if (ii == jj) {
            ind <- which(design_matrix[ii,] > 0.5)
            k2_grad[kk, ii] <- (2 * alpha_marg[ii] + 1) * kappa_marg[ii]
            k2_grad[kk, p + ind] <- alpha_marg[ii]^2 + alpha_marg[ii]
          } else {
            ind <- which(apply(design_matrix, 2, function(a) {
              b <- rep(0, p)
              b[c(ii, jj)] <- 1
              return (sum(abs(b - a)))
            }) < 0.5)
            if (length(ind) > 0.5) {
              k2_grad[kk, c(ii, jj)] <- theta_par[ind, 1] * c(alpha_marg[jj], alpha_marg[ii])
              k2_grad[kk, p + ind] <- alpha_marg[ii] * alpha_marg[jj]
            }
          }
          kk <- kk + 1
        }
      }
      grad <- 2 * colSums(matrix(k1_diff, nrow = p, ncol = k, byrow = FALSE) * k1_grad) +
        2 * colSums(matrix(k2_diff, nrow = p2, ncol = k, byrow = FALSE) * k2_grad)
      
      return (list("objective" = val, "gradient" = grad))
    }
    
    if (hasArg(x0)) x0 <- list(...)$x0 else x0 <- x0_negBin(k1_sample, k2_sample, design_matrix)
    optscontrol <- list(algorithm = "NLOPT_LD_MMA", xtol_rel = 1e-05, maxeval = 10000, 
                        print_level = 0, check_derivatives = FALSE)
    ub <- rep(Inf, k)
    sol <- nloptr::nloptr(x0 = x0, eval_f = obj, lb = rep(0, k), ub = ub, opts = optscontrol)
    levy_par <- theta2theta_par(sol$solution, design_matrix)
  }
  
  return (levy_par)
}

theta2theta_par <- function(theta, design_matrix) {
  # theta contains parameters (alpha_i..., m) - m is of length ncol(design_matrix)
  theta_par <- matrix(NA, nrow = ncol(design_matrix), ncol = 3)
  alpha <- theta[1:nrow(design_matrix)]
  theta_par[, 1] <- theta[-(1:nrow(design_matrix))]
  for (ii in 1:ncol(design_matrix)) {
    if (sum(design_matrix[, ii]) < 1.5) {
      ind <- which(design_matrix[, ii] > 0.5)
      theta_par[ii, 2] <- alpha[ind] / (1 + alpha[ind])
    } else {
      ind <- which(design_matrix[, ii] > 0.5)
      theta_par[ii, 2:3] <- alpha[ind] / (1 + alpha[ind])
    }
  }
  return (theta_par)
}

x0_negBin <- function(k1, k2, design_matrix) {
  p <- length(k1)
  ind <- rep(FALSE, length(k2))
  kk <- 1
  for (ii in 1:p) {
    for (jj in ii:p) {
      if (ii == jj) ind[kk] <- TRUE
      kk <- kk + 1
    }
  }
  k2 <- k2[ind]
  
  x0_mat <- matrix(NA, nrow = p, ncol = 2)
  for (ii in 1:p) {
    x0_mat[ii,] <- levy_cum_fit("negBin", k1[ii], k2[ii])
  }
  x0 <- x0_mat[, 2] / (1 - x0_mat[, 2])
  x0_m <- rep(0, ncol(design_matrix))
  x0_m[colSums(design_matrix) == 1] <- x0_mat[, 1]
  x0 <- c(x0, x0_m)
  return (x0)
}
