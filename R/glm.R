#' Loglikelihood adjustments for glm fits
#'
#' In a generalised linear model (glm), the user can choose between a range of
#' distributions of a response \eqn{y}, and can allow for non-linear relations
#' between the mean outcome for a particular combination of covariates \eqn{x},
#' \eqn{E(y_i\mid x_i)=\mu_i}, and the linear predictor,
#' \eqn{\eta_i=x_i^T\beta}, which is the link function \eqn{g(\mu_i)=\eta_i}.
#' it is required to be monotonic. (For a quick introduction, see
#' Kleiber and Zeileis (2008, Ch. 5.1), for more complete coverage of the
#' topic, see, for example, Davison (2003, Ch. 10.3))
#'
#' For more usage examples and more information on `glm` models, see the
#' *Introducing `chantrics`* vignette by running
#' \code{vignette("chantrics-vignette", package = "chantrics")}
#'
#' @section Supported families (within each family, any link function should work):
#'
#'   * `gaussian`
#'   * `poisson`
#'   * `binomial`
#'   * `MASS::negative.binomial`
#'
#'   Also works for [MASS::glm.nb()], note that the standard errors of the theta
#'   are not adjusted.
#'
#' @examples
#' # binomial example from Applied Econometrics in R, Kleiber/Zeileis (2008)
#' # ==  probit  ==
#' data("SwissLabor", package = "AER")
#' swiss_probit <- glm(participation ~ . + I(age^2),
#'   data = SwissLabor,
#'   family = binomial(link = "probit")
#' )
#' summary(swiss_probit)
#' swiss_probit_adj <- adj_loglik(swiss_probit)
#' summary(swiss_probit_adj)
#'
#' # == logit ==
#' swiss_logit <- glm(participation ~ . + I(age^2),
#'   data = SwissLabor,
#'   family = binomial(link = "logit")
#' )
#' summary(swiss_logit)
#' swiss_logit_adj <- adj_loglik(swiss_logit)
#' summary(swiss_logit_adj)
#' @references Davison, A. C. 2003. Statistical Models. Cambridge Series on
#'   Statistical and Probabilistic Mathematics 11. Cambridge University Press,
#'   Cambridge.
#'
#'   Kleiber, Christian, and Achim
#'   Zeileis. 2008. Applied Econometrics with R. Edited by Robert Gentleman,
#'   Kurt Hornik, and Giovanni Parmigiani. Use r! New York: Springer-Verlag.
#' @name glm
#' @aliases glm.nb
NULL

#' @export

# handling of dispersion parameters:
# http://people.stat.sfu.ca/~raltman/stat402/402L25.pdf

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
  response_vec <- get_response_from_model(object)
  df.resid <- df.residual(object)
  # theta for constant negbin models
  theta <- try(object$call$family$theta, silent = TRUE)

  if (inherits(object, "negbin")) {
    family <- "negbin"
    # get theta if saved
    if (!is.numeric(theta)) {
      theta <- rlang::env_get(bypasses.env, "theta", default = NULL)
      if (!is.numeric(theta)) {
        theta <- rlang::env_get(bypasses.env, paste0("theta", length(stats::coef(object))), default = NULL)
      }
    }
  } else {
    family <- object$family$family
  }
  llv <- glm_type_llv(
    family = family,
    x_mat = x_mat,
    pars = pars,
    response_vec = response_vec,
    linkinv = object$family$linkinv,
    df.resid = df.resid,
    theta = theta
  )

  # return other attributes from logLik objects
  if (inherits(object, "negbin")) {
    pars <- c(pars, "theta")
  }
  attr(llv, "df") <- length(pars)
  attr(llv, "nobs") <- stats::nobs(object)
  class(llv) <- "logLik_vec"
  return(llv)
}

glm_type_llv <- function(family, x_mat, pars, response_vec, linkinv, df.resid = NULL, theta = NULL, hurdle = c("no", "count", "zero")) {
  hurdle <- match.arg(hurdle)
  linkinv <- match.fun(linkinv)
  dimcheck(x_mat, pars)
  eta_vec <- x_mat %*% pars
  mu_vec <- linkinv(eta_vec)
  # "geometric" is negbin with theta = 1
  if (family == "geometric") {
    theta <- 1
    family <- "Negative Binomial("
  }
  if (family == "poisson") {
    llv <- stats::dpois(response_vec, lambda = mu_vec, log = TRUE)
  } else if (family == "binomial") {
    llv <- stats::dbinom(response_vec, 1, prob = mu_vec, log = TRUE)
  } else if (family == "gaussian") {
    # estimate the dispersion parameter
    disp <- dispersion.gauss(response_vec, mu_vec, df.resid - 1)
    pars <- c(pars, disp)
    names(pars)[length(pars)] <- "Dispersion"
    llv <- stats::dnorm(response_vec - mu_vec, mean = 0, sd = sqrt(disp), log = TRUE)
  } else if (substr(family, 1, 18) == "Negative Binomial(") {
    llv <- stats::dnbinom(response_vec, size = theta, mu = mu_vec, log = TRUE)
  } else if (family == "negbin") {
    # if theta should not be estimated (if through confint, etc.), then this
    # should be set previously.
    # If theta is not a number, reestimate
    if (!is.numeric(theta)) {
      theta <- MASS::theta.ml(y = response_vec, mu = mu_vec)
    }
    llv <- stats::dnbinom(response_vec, size = theta, mu = mu_vec, log = TRUE)
    rlang::env_poke(bypasses.env, "last_est_theta", theta)
  } else {
    rlang::abort(paste0(family, " is not supported."), class = "chantrics_not_supported_glm_family")
  }
  if (hurdle == "count") {
    # subtract the probability of the function being equal to 0,
    # see countregression-cam, Eqn. 4.54
    if (any(family == "negbin", family == "Negative Binomial(")) {
      llv <- llv - log(1 - stats::dnbinom(response_vec, size = theta, mu = mu_vec, log = FALSE))
    } else if (family == "poisson") {
      llv <- llv - log(1 - stats::dpois(response_vec, lambda = mu_vec, log = FALSE))
    } else {
      # safety break
      rlang::abort("Could not match count distribution specification.")
    }
    # set all contributions whose count is 0 to 0
    llv[response_vec == 0] <- 0
  }
  return(llv)
}

#' @keywords internal

dispersion.gauss <- function(response_vec, mu_vec, df) {
  # https://statmath.wu.ac.at/courses/heather_turner/glmCourse_001.pdf
  # (last slide in GLM > Estimation section)

  # http://people.stat.sfu.ca/~raltman/stat402/402L25.pdf (p. 4) -> dispersion
  # is the residual variance
  # final source: ISL, page 80, eqn. 3.25 (RSE formula)
  # response_vec are the true Y, the mu_vec are the fitted values

  return(sum((response_vec - mu_vec)^2) / (df))
}

#' @keywords internal

dispersion.stat <- function(response_vec, mu_vec, object) {
  # modelcount-hilbe, pg. 78 (Eqn. 3.4, pearson chi-sq statistic), pg. 79,
  # "The dispersion statistic of the Poisson model is defined as the
  # Pearson Chi2 statistic divided by the residual degrees of freedom"
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

#' @keywords internal

fittedhelper.glm <- function(object, type) {
  x_mat <- get_design_matrix_from_model(object)
  pars <- stats::coef(object)
  result <- x_mat %*% pars
  if (type == "response") {
    result <- get_unadj_object(object)$family$linkinv(result)
  }
  dim(result) <- NULL
  names(result) <- seq(1, length(result), 1)
  return(result)
}
