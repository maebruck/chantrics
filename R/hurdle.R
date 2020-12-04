#' Loglikelihood adjustments for pscl::hurdle fits
#'
#' Adjust the loglikelihood and the standard errors of a fitted [pscl::hurdle()] model.
#'
#' Describe zeroinfl models here
#'
#' Note that the [pscl::hurdle()] model has to be run with the option `x = TRUE` in order for the adjustment to execute properly
#'
#' @section Supported families (within each family, any link function should work):
#'
#' * `geometric`
#' * `poisson`
#' * `negbin`
#' * `binomial` (for the zero mass distribution only)
#'
#' @examples
#' # hurdle model from AER, pg. 139-140
#' library(pscl)
#' data("RecreationDemand", package = "AER")
#' rd_hurdle <- hurdle(trips ~ . | quality + income, data = RecreationDemand, dist = "negbin", x = TRUE)
#' summary(rd_hurdle)
#'
#' # adjust model
#' adj_loglik(rd_hurdle)
#'
#'
#' @name hurdle
NULL

#' @export


logLik_vec.hurdle <- function(object, pars = NULL, ...) {
  if (!missing(...)) {
    rlang::warn("extra arguments discarded")
  }
  # import coefficients
  if (is.null(pars)) {
    pars <- stats::coef(object)
  }
  # issue warning if the model has not converged
  if (!object$converged) {
    rlang::warn("The original model's parameters did not converge.")
  }
  # calculate mus
  # try getting the design matrix from object
  count_mat <- get_design_matrix_from_model(object, "count")
  zero_mat <- get_design_matrix_from_model(object, "zero")
  #split the parameters into the two models
  count_pars <- pars[startsWith(names(pars), "count")]
  zero_pars <- pars[startsWith(names(pars), "zero")]
  response_vec <- get_response_from_model(object)


  # return other attributes from logLik objects
  attr(llv, "df") <- length(pars)
  attr(llv, "nobs") <- stats::nobs(object)
  class(llv) <- "logLik_vec"
  return(llv)
}

#' @keywords internal

dispersion.gauss <- function(response_vec, mu_vec, df) {
  # https://statmath.wu.ac.at/courses/heather_turner/glmCourse_001.pdf (last slide in GLM > Estimation section)

  # http://people.stat.sfu.ca/~raltman/stat402/402L25.pdf (p. 4) -> dispersion is the residual variance
  # final source: ISL, page 80, eqn. 3.25 (RSE formula)
  # response_vec are the true Y, the mu_vec are the fitted values

  return(sum((response_vec - mu_vec)^2) / (df))
}

#' @keywords internal

dispersion.stat <- function(response_vec, mu_vec, object) {
  # modelcount-hilbe, pg. 78 (Eqn. 3.4, pearson chi-sq statistic), pg. 79, "The dispersion statistic of the Poisson model is defined as the Pearson Chi2 statistic divided by the residual degrees of freedom"
  return(sum(((response_vec - mu_vec)^2) / object$family$variance(mu_vec)) / (stats::df.residual(object)))
}

#' @keywords internal

dimcheck <- function(x_mat, pars) {
  # check that dimensions of xmat and pars coincide
  if (ncol(x_mat) != length(pars)) {
    rlang::abort(paste0(
      "The length of 'pars' (", length(pars), ") does not fit\n",
      "the number of columns in the design matrix, ", ncol(x_mat)
    ),
    class = "chantrics_pars_wrong_length"
    )
  }
  if (any(colnames(x_mat) != names(pars))) {
    rlang::warn(paste0(
      "names(pars) is not equal to the colnames() of the design matrix.\n",
      "Continuing with the parameters as matched below:",
      "Design matrix: ", paste(colnames(x_mat), sep = ", "),
      "\n names(pars): ", paste(names(pars), sep = ", ")
    ),
    class = "chantrics_parnames_do_not_match"
    )
  }
}

#' @export

logLik_vec.negbin <- function(object, pars = NULL, ...) {
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
  dimcheck(x_mat, pars)
  eta_vec <- x_mat %*% pars
  mu_vec <- object$family$linkinv(eta_vec)
  # try getting the response vector from glm
  response_vec <- get_response_from_model(object)
  # theta <- dispersion.stat(response_vec, mu_vec, object)
  # bypass theta calculation
  theta <- try(get("theta", envir = bypasses.env), silent = TRUE)
  if (is.error(theta)) {
    theta <- MASS::theta.ml(y = response_vec, mu = mu_vec)
  }
  pars <- c(pars, "theta")
  llv <- stats::dnbinom(response_vec, size = theta, mu = mu_vec, log = TRUE)

  # return other attributes from logLik objects
  attr(llv, "df") <- length(pars)
  attr(llv, "nobs") <- stats::nobs(object)
  class(llv) <- "logLik_vec"
  return(llv)
}

#' @keywords internal

fittedhelper.glm <- function(object, modelname) {
  x_mat <- get_design_matrix_from_model(object)
  pars <- stats::coef(object)
  eta_vec <- x_mat %*% pars
  mu_vec <- get_unadj_object(object)$family$linkinv(eta_vec)
  dim(mu_vec) <- NULL
  names(mu_vec) <- seq(1, length(mu_vec), 1)
  return(mu_vec)
}
