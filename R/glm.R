#' Loglikelihood adjustments for glm fits
#'
#' Describe glm methods here
#'
#' @section Supported families (within each family, any link function should work):
#'
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
  x_mat <- try(stats::model.matrix(object), silent = TRUE)
  if (is.error(x_mat)) {
    x_mat <- try(object$x, silent = TRUE)
    if (is.error(x_mat)) {
      raise_yield_error(
        "glm object",
        "design matrix",
        "executing the glm with the option 'y = TRUE"
      )
    }
  }
  eta_vec <- x_mat %*% pars
  mu_vec <- object$family$linkinv(eta_vec)
  # try getting the response vector from glm
  response_vec <- try(stats::model.response(object), silent = TRUE)
  if (is.error(response_vec)) {
    response_vec <- try(object$y, silent = TRUE)
    if (is.error(response_vec)) {
      raise_yield_error(
        "glm object",
        "response vector",
        "executing the glm with the option 'y = TRUE'"
      )
    }
  }
  if (object$family$family == "poisson") {
    llv <- stats::dpois(response_vec, lambda = mu_vec, log = TRUE)
  } else if (object$family$family == "binomial") {
    llv <- stats::dbinom(response_vec, 1, prob = mu_vec, log = TRUE)
  }
  # add other densities here
  # return other attributes from logLik objects
  attr(llv, "df") <- length(pars)
  attr(llv, "nobs") <- stats::nobs(object)
  class(llv) <- "logLik_vec"
  return(llv)
}
