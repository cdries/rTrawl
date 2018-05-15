#' Estimation of Trawl processes
#'
#' display variance signature plot
#'
#' TODO
#'
#' CITE TODO.
#' @name vs_trawl
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

  }



  return(vs)
}
