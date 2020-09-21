#' Evaluate loglikelihood contributions from specific observations
#'
#' Generic function for calculating loglikelihood contributions from
#' individual observations for a fitted model.
#'
#' @param object A fitted model object.
#' @param ... Further arguments.
#' @export

logLik_vec <- function(object, ...) {
  UseMethod("logLik_vec")
}

#' Log-likelihood adjustments for fitted models
#'
#' This function adjusts the log-likelihood of fitted model objects based on
#' \href{http://doi.org/10.1093/biomet/asm015}{Chandler and Bate (2007)}.
#'
#' @export

# if required, turn this into a method (see logLik_vec) and the below into the
# .default() method

adj_loglik <- function(x,
                       cluster = NULL,
                       use_vcov = TRUE,
                       ...) {
  #check if x is a supported model type
  #adjust x
  adjusted_x <-
    chant_obj(x, cluster = cluster, use_vcov = use_vcov, ...)
  class(adjusted_x) <- c("chantrics", "chandwich", class(x))
  return(adjusted_x)
}
