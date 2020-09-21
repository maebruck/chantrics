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
#' ANOVA tables: compare nested models
#'
#' \code{anova} method for \code{chantrics} objects
#'
#' Compare two or more nested models that have been adjusted using the
#' \code{\link{adj_logLik}} method. It uses the adjusted likelihood ratio test
#' statistic (ALRTS), as described in Section 3.5 of
#' \href{http://doi.org/10.1093/biomet/asm015}{Chandler and Bate (2007)}.
#'
#' Note that the objects specified do not have to be sorted in a specific way,
#' they will be sorted by the function by the number of parameters, as returned
#' by \code{attr(model1, "p_current")}.
#'
#' @param model1 Object of class \code{chantrics}, as returned by
#'   \code{\link{adj_logLik}}.
#' @param model2 Object of class \code{chantrics}, as returned by
#'   \code{\link{adj_logLik}}.
#' @param ... Further objects of class \code{chantrics}, as returned by
#'   \code{\link{adj_logLik}}, and/or parameters that will be passed to
#'   \code{\link[chandwich]{anova.chandwich}}, and then further to
#'   \code{\link[chandwich]{compare_models}}. The type of adjustment, out of
#'   \code{"vertical"}, \code{"cholesky"}, \code{"spectral"}, \code{"none"}, as
#'   specified in the parameter \code{type}, can also be specified here.
#'
#' @export
anova.chantrics <- function(model1, model2, ...) {
  dotargs <- list(...)
  potential_model_objects <-
    c(model1, model2, subset(dotargs, names(dotargs) == ""), use.names = FALSE)
  #check whether the unnamed objects are actually chantrics objects
  pmo_is_chantrics <- sapply(potential_model_objects, function(x){is(x, "chantrics")})
  #warn user that we drop supplied unnamed arguments that are not chantrics objects
  if(any(pmo_is_chantrics)){
    warn("One or more of the unnamed objects supplied are not 'chantrics' objects,\nwhich have been dropped.")
  }
  #drop non-chantrics objects
  model_objects <- subset(potential_model_objects, pmo_is_chantrics)
  model_objects
  #find out if sorting is necessary
  #prepare for piping into chandwich::anova.chandwich
}
