#' Evaluate loglikelihood contributions from specific observations
#'
#' Generic function for calculating the loglikelihood contributions from
#' individual observations for a fitted model.
#'
#' @param object A fitted model object.
#' @param ... Further arguments.
#'
#' @return An object of class `"logLik_vec"`, which is a numeric vector of length
#'   `nobs(object)` (i.e. the number of observations in `object`) of the
#'   loglikelihood of each observation. Additionally, it contains the
#'   attributes `df` (model degrees of freedom) and `nobs` (number of
#'   observations).
#'
#'   The methods [stats::logLik()], and [stats::nobs()] are available.
#'
#' @seealso [stats::logLik()]
#'
#' @export

logLik_vec <- function(object, ...) {
  UseMethod("logLik_vec")
}

#' @importFrom stats logLik
#' @export

logLik.logLik_vec <- function(object, ...) {
  if (!missing(...)) {
    rlang::warn("extra arguments discarded")
  }
  val <- sum(object)
  # logLik_vec has "nobs" and "df" attributes -> just pass through
  attributes(val) <- attributes(object)[c("nobs", "df")]
  class(val) <- "logLik"
  return(val)
}

#' @importFrom stats nobs
#' @export


nobs.logLik_vec <- function(object, ...) {
  if (!missing(...)) {
    rlang::warn("extra arguments discarded")
  }
  return(attr(object, "nobs"))
}
