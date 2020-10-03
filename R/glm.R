#' Loglikelihood adjustments for glm fits
#'
#' Describe glm methods here
#'
#' @section Supported families (within each family, any link function should work):
#'
#' * `poisson`
#'
#' @name glm
NULL

#' @export

logLik_vec.glm <- function(object, pars = NULL, ...) {
  if (!missing(...)) {
    rlang::warn("extra arguments discarded")
  }
  #import coefficients
  if (is.null(pars)) {
    pars <- stats::coef(object)
  }
  #calculate mus
  #try getting the design matrix from glm object
  x_mat <- try(stats::model.matrix(object), silent = TRUE)
  if (is.error(x_mat)) {
    x_mat <- try(object$x, silent = TRUE)
    if (is.error(x_mat)) {
      raise_yield_error("glm object",
                        "design matrix",
                        "executing the glm with the option 'y = TRUE")
    }
  }
  eta_vec <- x_mat %*% pars
  mu_vec <- object$family$linkinv(eta_vec)
  #try getting the response vector from glm
  response_vec <- try(stats::model.response(object), silent = TRUE)
  if (is.error(response_vec)) {
    response_vec <- try(object$y, silent = TRUE)
    if (is.error(response_vec)) {
      raise_yield_error("glm object",
                        "response vector",
                        "executing the glm with the option 'y = TRUE'")
    }
  }
  if (object$family$family == "poisson") {
    llv <- stats::dpois(response_vec, lambda = mu_vec, log = TRUE)
  }
  #add other densities here
  #return other attributes from logLik objects
  attr(llv, "df") <- length(pars)
  attr(llv, "nobs") <- stats::nobs(object)
  class(llv) <- "logLik_vec"
  return(llv)
}
