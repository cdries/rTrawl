
bootstrap_trawl <- function(object, n_bootstrap = 1000, level = 0.95, univariate = TRUE, ...) {
  
  # first draw and fit
  sim <- sim_trawl(object, univariate = univariate, ...)
  
  sim$h <- object$h
  sim$include_b <- object$include_b
  sim$include_cum1 <- object$include_cum1
  sim$method <- object$method
  ft <- fit_trawl(sim)
  
  # make objects to save bootstrap results
  if (univariate) {
    trawl_est_all <- matrix(NA, nrow = n_bootstrap, ncol = length(ft$trawl_par))
    trawl_est_all[1,] <- ft$trawl_par
  } else {
    # TODO
  }
  
  for (ii in 2:n_bootstrap) {
    sim <- sim_trawl(object, univariate = univariate, ...)
    
    sim$h <- object$h
    sim$include_b <- object$include_b
    sim$include_cum1 <- object$include_cum1
    sim$method <- object$method
    ft <- fit_trawl(sim)
    
    if (univariate) trawl_est_all[ii,] <- ft$trawl_par # TODO else
  }
  
  # construct confidence intervals
  a2 <- 0.5 * (1 - level)
  if (univariate) {
    ci <- apply(trawl_est_all, 2, function(a) sort(a)[c(max(1, floor(a2 * n_bootstrap)), 
                                                        min(n_bootstrap, ceiling((1 - a2) * n_bootstrap)))])
    summ <- cbind(object$trawl_par, t(ci))
    colnames(summ) <- c("estimate", "lb", "ub")
  } else {
    # TODO
  }
  
  return (list("summary" = summary, "estimates" = trawl_est_all))
}
