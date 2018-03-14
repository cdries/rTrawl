#' Sample ACF of irregularly observed time series
#'
#' sample autocorrelation of a trawl process
#'
#' TODO
#'
#' CITE TODO.
#' @name acf_sample
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
#' @export acf_sample
#' @useDynLib rTrawl
acf_sample <- function(object, h, dff = 0, lag_max = 25, drop_zero = TRUE, multi = 5) {
  
  x_grid <- object$x_grid
  p_grid <- object$p_grid
  T0 <- object$T0
  TT <- object$TT
  
  if (dff < 0.5) {
    # ACF of the process itself
    acfh <- as.numeric(acf_sample_p(h, x_grid, p_grid, T0, TT, lag_max, multi)) # TODO
  } else {
    # ACF of the differenced process
    acfh <- as.numeric(acf_sample_dp(h, x_grid, p_grid, T0, TT, lag_max, multi)) # TODO
  }
  
  if (!drop_zero) acft <- c(1, acft)
  
  return (acft)
  
  
  
  # T0_offset <- seq(0, h, length.out = multi + 1)[-(multi + 1)]
  # 
  # acft <- rep(0, length(diff_vec))
  # rv <- 0
  # 
  # for (mm in 1:multi) {
  #   obs_process <- observe_process(x_grid, p_grid, T0 + T0_offset[mm], TT, h)
  #   
  #   if (dff < 0.5) {
  #     obs_x <- as.integer(round((obs_process$x_grid - T0 - T0_offset[mm]) / h))
  #     obs_p <- obs_process$p_grid
  #     n <- length(obs_p)
  #     n_total <- floor((TT - T0 - T0_offset[mm]) / h)
  #     obs_p <- obs_p - (obs_p[n] * (n_total - obs_x[n] + 1) + sum(obs_p[-n] * diff(obs_x))) / n_total
  #     for (ii in 1:length(acft)) {
  #       acft[ii] <- acft[ii] +
  #         .Call('acf_helper', as.integer(c(obs_x, n_total)), c(obs_p, obs_p[n]),
  #               length(obs_x) + 1, diff_vec[ii], PACKAGE = "rTrawl")
  #     }
  #     rv <- rv + (obs_p[n]^2 * (n_total - obs_x[n] + 1) + sum(obs_p[-n]^2 * diff(obs_x)))
  #     
  #   } else {
  #     
  #     obs_x <- as.integer(round((obs_process$x_grid[-1] - T0 - T0_offset[mm]) / h))
  #     obs_dp <- diff(obs_process$p_grid)
  #     obs_dp <- obs_dp - sum(obs_dp) / floor((TT - T0 - T0_offset[mm]) / h - 1)
  #     for (ii in 1:length(acft)) {
  #       ind_current <- (obs_x - diff_vec[ii]) %in% obs_x
  #       ind_lagged <- (obs_x + diff_vec[ii]) %in% obs_x
  #       acft[ii] <- acft[ii] + sum(obs_dp[ind_current] * obs_dp[ind_lagged]) / floor((TT - T0) / h - 1 - diff_vec[ii])
  #     }
  #     rv <- rv + sum(obs_dp^2) / floor((TT - T0 - T0_offset[mm]) / h - 1)
  #   }
  # }
  # acft <- acft / rv
  # 
  # if (!drop_zero) acft <- c(1, acft)
  
  return (acft)
}


#' ACF of Trawl processes
#'
#' theoretical autocorrelation of a trawl process
#'
#' TODO
#'
#' CITE TODO.
#' @name acf_trawl
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
#' @export acf_trawl
acf_trawl <- function(object, h, dff = 0, lag_max = 25, drop_zero = TRUE) {
  
  # h_vec <- (0:lag_max) * h
  trawl <- object$trawl
  trawl_par <- object$trawl_par
  b <- object$b
  
  if (dff < 0.5) {
    # ACF of the process itself - 
    # TODO: make correct when b not equal to zero
    acfh <- as.numeric(acf_trawl_p(h, trawl, trawl_par, lag_max))
  } else {
    # ACF of the differenced process
    acfh <- as.numeric(acf_trawl_dp(h, trawl, trawl_par, b, lag_max))
  }
  
  if (!drop_zero) acft <- c(1, acft)
  
  return (acft)
}
