#' Loglikelihood adjustments for GLM fits
#'
#' Describe glm methods here
#'
#' @name glm
NULL

#' @rdname evd
#' @export

alogLik.glm <- function(x, cluster = NULL, use_vcov = NULL, ...){
  #forward into specific functions
  x$method <- glm_chand_fit
  if (x$family$family == "poisson" && x$family$link == "log"){
    alogLik.glm.poisson(x, cluster = cluster, use_vcov = use_vcov, ...)
  } else {
    stop("This combination of link function and family is currently not supported.")
  }
}

#' Custom function for glm to fit the model using the adjustment
#'
#' @export

glm_chand_fit <- function(x, y, ){
  adjust_loglik
}

#' @export

logLikVec.glm <- function(fitted_object, pars) {
  #calculate mus
  #try getting the design matrix from glm object
  x_mat <- try(model.matrix(fitted_object),silent = TRUE)
  if (is.error(x_mat)){
    x_mat <- try(fitted_object$x, silent = TRUE)
    if (is.error(x_mat)){
      raise_yield_error("glm object", "design matrix", "executing the glm with the option 'y = TRUE")
    }
  }
  eta_vec <- x_mat %*% pars
  mu_vec <- fitted_object$family$linkinv(eta_vec)
  #try getting the response vector from glm
  response_vec <- try(model.response(fitted_object), silent = TRUE)
  if (is.error(response_vec)){
    response_vec <- try(fitted_object$y, silent = TRUE)
    if (is.error(response_vec)){
      raise_yield_error("glm object", "response vector", "executing the glm with the option 'y = TRUE'")
    }
  }
  if (fitted_object$family$family == "poisson"){
    return(dpois(response_vec, lambda = mu_vec, log = TRUE))
  }
}


