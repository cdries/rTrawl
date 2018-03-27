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

  # trawl measure
  trawl <- object$trawl
  trawl_par <- object$trawl_par
  
  # for multivariate fit
  design_matrix <- object$design_matrix
  
  # cumulants - compute if required
  if ("k1_sample" %in% names(object)) {
    k1_sample <- object$k1_sample
    k2_sample <- object$k2_sample
  } else {
    T0 <- object$T0
    TT <- object$TT
    x_grid <- object$x_grid
    p_grid <- object$p_grid
    if (length(trawl) == 1) {
      lebA <- leb_AtA(0.0, trawl, trawl_par)
      k1_sample <- cum_sample(1L, as.numeric(x_grid), as.numeric(p_grid), as.numeric(TT)) / lebA
      k2_sample <- cum_sample(2L, as.numeric(x_grid), as.numeric(p_grid), as.numeric(TT)) / lebA
    } else {
      p <- length(trawl)
      k1_sample <- rep(NA, p)
      k2_sample <- rep(NA, p * (p + 1) / 2)
      kk <- 1
      for (ii in 1:p) {
        lebAii <- leb_AtA(0.0, trawl[[ii]], trawl_par[[ii]])
        k1_sample[ii] <- cum_sample(1L, as.numeric(x_grid[[ii]]), as.numeric(p_grid[[ii]]), as.numeric(TT)) / lebAii
        for (jj in ii:p) {
          lebAjj <- leb_AtA(0.0, trawl[[jj]], trawl_par[[jj]])
          k2_sample[kk] <- ccf_sample_p(0.0, x_grid[[ii]], p_grid[[ii]], x_grid[[jj]], p_grid[[jj]], TT, 0L) / sqrt(lebAii * lebAjj)
          # TODO: scale to correct mqgnitude!
          kk <- kk + 1
        }
      }
    }
  }
  
  # fit using 1st (and 2nd) cumulant
  if (NROW(design_matrix) == 1) {
    lfit <- levy_cum_fit(levy_seed, k1_sample, k2_sample)
  } else {
    constraints <- diag(ncol(design_matrix))
    lfit <- levy_cum_fit_mv(levy_seed, k1_sample, k2_sample, design_matrix, constraints)
  }
  
  return(lfit)
}


levy_cum_fit_mv <- function(levy_seed, k1_sample, k2_sample, design_matrix, constraints) {
  
  p <- nrow(design_matrix)
  p2 <- p * (p + 1) / 2
  k <- ncol(constraints)

  if (levy_seed == "Poisson" || levy_seed == "Skellam") {
    
    obj <- function(theta) {
      
      theta <- constraints %*% theta
      
      k1_theor <- design_matrix %*% theta
      k2_theor <- (design_matrix %*% diag(theta) %*% t(design_matrix))[lower.tri(diag(p), diag = TRUE)]
      
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
    
    optscontrol <- list(algorithm = "NLOPT_LD_MMA", xtol_rel = 1e-05, maxeval = 10000, 
                        print_level = 3, check_derivatives = FALSE)
    sol <- nloptr::nloptr(x0 = x0, eval_f = obj, lb = rep(0, k), ub = rep(Inf, k), opts = optscontrol)
    levy_par <- matrix(constraints %*% sol$solution, ncol = 1)
  }
  
  return (levy_par)
}


