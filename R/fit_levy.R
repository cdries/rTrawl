#' Estimation of Trawl processes
#'
#' estimates a path of a Trawl process
#'
#' TODO
#'
#' CITE TODO.
#' @name fit_levy
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
      k1_sample <- cum_sample(1L, as.numeric(x_grid), as.numeric(p_grid), as.numeric(TT)) / lebA
      k2_sample <- cum_sample(2L, as.numeric(x_grid), as.numeric(p_grid), as.numeric(TT)) / lebA
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
    
    if (hasArg(x0)) x0 <- list(...)$x0 else x0 <- rep(0.1, k)
    optscontrol <- list(algorithm = "NLOPT_LD_MMA", xtol_rel = 1e-05, maxeval = 10000, 
                        print_level = 0, check_derivatives = FALSE)
    ub <- rep(Inf, k)
    # ub[1:nrow(design_matrix)] <- 1
    sol <- nloptr::nloptr(x0 = x0, eval_f = obj, lb = rep(0, k), ub = ub, opts = optscontrol)
    levy_par <- theta2theta_par(sol$solution, design_matrix)
  }
  
  return (levy_par)
}

# theta_par2theta_marg <- function(levy_par, design_matrix) {
#   theta_marg <- rep(0, nrow(design_matrix))
#   for (ii in 1:ncol(design_matrix)) {
#     if (sum(design_matrix[, ii]) < 1.5) {
#       ind <- which(design_matrix[, ii] > 0.5)
#       theta_marg[ind] <- theta_marg[ind] + levy_par[ii, 2]
#     } else {
#       ind <- which(design_matrix[, ii] > 0.5)
#       theta_marg[ind] <- theta_marg[ind] + levy_par[ii, 2:3]
#     }
#   }
#   theta_marg <- theta_marg / rowSums(design_matrix)
#   return (theta_marg)
# }

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
