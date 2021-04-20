#' Loglikelihood adjustments for pscl::hurdle fits
#'
#' NOTE: *Hurdle models are currently not supported.*
#'
#' Adjust the loglikelihood and the standard errors of a fitted [pscl::hurdle()]
#' model.
#'
#' Note that the [pscl::hurdle()] model has to be run with the option `x = TRUE`
#' in order for the adjustment to execute properly. The functions
#' [residuals.chantrics()] and [`fitted.chantrics()`][fitted()] are currently
#' disabled for `hurdle` models. Additionally, sequential [anova.chantrics()]
#' are not available.
#'
#' @section Supported families:
#' Within each family, any link function should work.
#'
#'   * `geometric`
#'   * `poisson`
#'   * `negbin`
#'   * `binomial` (for the zero mass distribution only)
#'
#'
#' @examples
#' # hurdle model from AER, pg. 139-140
#' library(pscl)
#' data("RecreationDemand", package = "AER")
#' rd_hurdle <- hurdle(trips ~ . | quality + income,
#'   data = RecreationDemand,
#'   dist = "negbin", x = TRUE
#' )
#' summary(rd_hurdle)
#'
#' # adjust model
#' # adj_loglik(rd_hurdle)
#' @name hurdle
#' @keywords internal
NULL

#' @export


logLik_vec.hurdle <- function(object, pars = NULL, ...) {
  if (!missing(...)) {
    rlang::warn("extra arguments discarded")
  }
  # import coefficients
  if (any(is.null(pars), length(pars) == 0)) {
    pars <- stats::coef(object)
  }
  # issue warning if the model has not converged
  if (!object$converged) {
    rlang::warn("The original model's parameters did not converge.")
  }
  # if none is saved, reestimate
  theta <- rlang::env_get(bypasses.env, "theta", default = NULL)
  if (is.null(theta)) {
    count_theta <- NULL
    zero_theta <- NULL
  } else {
    # if not error, parse
    if (!is.null(theta[["count"]])) {
      count_theta <- theta[["count"]]
    } else {
      count_theta <- NULL
    }
    if (!is.null(theta[["zero"]])) {
      zero_theta <- theta[["zero"]]
    } else {
      zero_theta <- NULL
    }
  }

  response_vec <- get_response_from_model(object)
  # calculate mus
  # try getting the design matrix from object
  count_mat <- get_design_matrix_from_model(object, "count")
  # add "count_" to the colnames in order for pars to match
  colnames(count_mat) <- vapply(
    colnames(count_mat),
    function(x) paste0("count_", x),
    character(1L)
  )
  # split the parameters into the two models
  count_pars <- try(pars[startsWith(names(pars), "count")], silent = TRUE)
  if (is.error(count_pars)) {
    # if pars has no names, reconstruct from object$coefficients
    count_mle_names <- names(object[["coefficients"]][["count"]])
    count_pars <- pars[seq_along(count_mle_names)]
    names(count_pars) <- vapply(
      count_mle_names,
      function(x) paste0("count_", x),
      character(1L)
    )
  }
  count_family <- object$dist$count
  count_linkinv <- function(eta) {
    pmax(exp(eta), .Machine$double.eps)
  }
  count_df.resid <- length(count_pars)
  if (count_family == "negbin") {
    count_df.resid <- count_df.resid + 1
  }
  # remove all observations that are 0
  count_llv <- glm_type_llv(
    family = count_family,
    x_mat = count_mat,
    pars = count_pars,
    response_vec = response_vec,
    linkinv = count_linkinv,
    df.resid = count_df.resid,
    theta = count_theta,
    hurdle = "count"
  )
  # print(sum(count_llv))
  count_theta_est <- rlang::env_get(bypasses.env, "last_est_theta",
    default = NULL
  )
  if (rlang::env_has(bypasses.env, "last_est_theta")) {
    rlang::env_unbind(bypasses.env, "last_est_theta")
  }


  zero_mat <- get_design_matrix_from_model(object, "zero")
  # add "count_" to the colnames in order for pars to match
  colnames(zero_mat) <- vapply(
    colnames(zero_mat),
    function(x) paste0("zero_", x),
    character(1L)
  )
  # split the parameters into the two models
  zero_pars <- try(pars[startsWith(names(pars), "zero")])
  if (is.error(zero_pars)) {
    # if pars has no names, reconstruct from object$coefficients
    zero_mle_names <- names(object[["coefficients"]][["zero"]])
    zero_pars <- pars[(length(pars) - length(zero_mle_names)):length(pars)]
    names(zero_pars) <- vapply(
      zero_mle_names,
      function(x) paste0("zero_", x),
      character(1L)
    )
  }
  zero_family <- object$dist$zero
  zero_linkinv <- object$linkinv
  if (is.null(zero_linkinv)) {
    zero_linkinv <- function(eta) {
      pmax(exp(eta), .Machine$double.eps)
    }
  }
  zero_df.resid <- length(zero_pars)
  # zero_theta <- try(object$theta[["zero"]], silent = TRUE)
  zero_response <- response_vec
  zero_response[zero_response != 0] <- 1
  zero_theta <- NULL
  if (zero_family == "binomial") {
    zero_llv <- glm_type_llv(
      family = zero_family,
      x_mat = zero_mat,
      pars = zero_pars,
      response_vec = zero_response,
      linkinv = zero_linkinv,
      df.resid = zero_df.resid,
      theta = zero_theta,
      hurdle = "zero"
    )
  } else {
    if (zero_family == "geometric") {
      zero_theta <- 1
    } else {
      zero_df.resid <- zero_df.resid + 1
    }
    null_eta_vec <- zero_mat %*% zero_pars
    null_mu_vec <- zero_linkinv(null_eta_vec)
    if (!is.numeric(zero_theta)) {
      zero_theta <- MASS::theta.ml(y = zero_response, mu = null_mu_vec)
    }
    lv_at_zero <- stats::dnbinom(rep(0, length(zero_response)),
      size = zero_theta,
      mu = null_mu_vec,
      log = FALSE
    )
    zero_llv <- (1 - zero_response) * log(lv_at_zero) + zero_response * log(1 - lv_at_zero)
  }
  # print(sum(zero_llv))
  llv <- count_llv + zero_llv
  # save theta objects in bypasses
  if (any("negbin" == c(count_family, zero_family))) {
    theta_ret <- list()
    if (is.numeric(count_theta_est)) {
      theta_ret[["count"]] <- count_theta_est
    }
    if (is.numeric(zero_theta)) {
      theta_ret[["zero"]] <- zero_theta
    }
    rlang::env_poke(bypasses.env, "negbin_theta_est", theta_ret)
  }

  # return other attributes from logLik objects

  attr(llv, "df") <- count_df.resid + zero_df.resid
  attr(llv, "nobs") <- stats::nobs(object)
  class(llv) <- "logLik_vec"
  return(llv)
}

#' @importFrom stats nobs
#' @export

nobs.hurdle <- function(object, ...) {
  return(object$n)
}
