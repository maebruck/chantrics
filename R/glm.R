#' Loglikelihood adjustments for glm fits
#'
#' Describe glm methods here
#'
#' @section Supported families (within each family, any link function should work):
#'
#' * `gaussian`
#' * `poisson`
#' * `binomial`
#'
#' @examples
#' # binomial example from Applied Econometrics in R, Kleiber/Zeileis (2008)
#' # ==  probit  ==
#' data("SwissLabor", package = "AER")
#' swiss_probit <- glm(participation ~ . + I(age^2), data = SwissLabor,
#'                     family = binomial(link = "probit"))
#' summary(swiss_probit)
#' swiss_probit_adj <- adj_loglik(swiss_probit)
#' summary(swiss_probit_adj)
#'
#' # == logit ==
#' swiss_logit <- glm(participation ~ . + I(age^2), data = SwissLabor,
#'                    family = binomial(link = "logit"))
#' summary(swiss_logit)
#' swiss_logit_adj <- adj_loglik(swiss_logit)
#' summary(swiss_logit_adj)
#' @name glm
NULL

#' @export

# handling of dispersion parameters: http://people.stat.sfu.ca/~raltman/stat402/402L25.pdf

logLik_vec.glm <- function(object, pars = NULL, ...) {
  if (!missing(...)) {
    rlang::warn("extra arguments discarded")
  }
  # import coefficients
  if (is.null(pars)) {
    pars <- stats::coef(object)
  }
  # calculate mus
  # try getting the design matrix from glm object
  x_mat <- get_design_matrix_from_model(object)
  eta_vec <- x_mat %*% pars
  mu_vec <- object$family$linkinv(eta_vec)
  # try getting the response vector from glm
  response_vec <- get_response_from_model(object)
  if (object$family$family == "poisson") {
    llv <- stats::dpois(response_vec, lambda = mu_vec, log = TRUE)
  } else if (object$family$family == "binomial") {
    llv <- stats::dbinom(response_vec, 1, prob = mu_vec, log = TRUE)
  } else if (object$family$family == "gaussian") {
    #estimate the dispersion parameter
    disp <- dispersion.gauss(response_vec, mu_vec, df.residual(object) - 1)
    pars <- c(pars, disp)
    names(pars)[length(pars)] <- "Dispersion"
    llv <- stats::dnorm(response_vec - mu_vec, mean = 0, sd = sqrt(disp), log = TRUE)
  }
  # add other densities here

  # return other attributes from logLik objects
  attr(llv, "df") <- length(pars)
  attr(llv, "nobs") <- stats::nobs(object)
  class(llv) <- "logLik_vec"
  return(llv)
}

#' @keywords internal

dispersion.gauss <- function(response_vec, mu_vec, df) {
  #https://statmath.wu.ac.at/courses/heather_turner/glmCourse_001.pdf (last slide in GLM > Estimation section)

  #http://people.stat.sfu.ca/~raltman/stat402/402L25.pdf (p. 4) -> dispersion is the residual variance
  #final source: ISL, page 80, eqn. 3.25 (RSE formula)
  #response_vec are the true Y, the mu_vec are the fitted values

  return(sum((response_vec - mu_vec) ^ 2) / (df))
}
