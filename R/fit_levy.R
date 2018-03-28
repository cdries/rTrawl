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
    if ("constraints" %in% names(object)) constraints <- object$constraints else constraints <- diag(ncol(design_matrix))
    lfit <- levy_cum_fit_mv(levy_seed, k1_sample, k2_sample, design_matrix, constraints, ...)
  }
  
  return(lfit)
}


levy_cum_fit_mv <- function(levy_seed, k1_sample, k2_sample, design_matrix, constraints, ...) {
  
  p <- nrow(design_matrix)
  p2 <- p * (p + 1) / 2
  k <- ncol(constraints)

  if (levy_seed == "Poisson" || levy_seed == "Skellam") {
    
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
  }
  
  return (levy_par)
}


