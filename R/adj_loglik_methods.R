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
