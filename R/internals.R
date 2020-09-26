#' chandwich internal functions
#'
#' Internal functions: These are not to be exported.
#' @name internal
#' @keywords internal
NULL

#' @keywords internal
#' @rdname internal

chant_obj <- function(x, cluster = NULL, ...) {
  #add tests if object valid
  #does x have a logLik_vec method?
  if (!any(paste0("logLik_vec.", class(x)) %in% utils::methods("logLik_vec"))) {
    rlang::abort("x does not have a logLik_vec method")
  }
  #create function for log-likelihood of x
  logLik_f <- function(pars, fitted_object, ...) {
    return(c(logLik_vec(fitted_object, pars = pars)))
  }
  name_pieces <- c(class(x))
  #add glm family to name
  try(name_pieces <- c(x$family$family, name_pieces), silent = TRUE)
  #get mle estimate from x
  mle = stats::coef(x)
  #adjust object using chandwich
  adjusted_obj <-
    chandwich::adjust_loglik(
      loglik = logLik_f,
      cluster = cluster,
      fitted_object = x,
      p = length(mle),
      par_names = names(mle),
      name = paste(name_pieces, collapse = "_"),
      mle = mle
    )
  class(adjusted_obj) <- c("chantrics", "chandwich")
  return(adjusted_obj)
}

#' @keywords internal
#' @rdname internal
is.error <- function(x)
  inherits(x, "try-error")

#' @keywords internal
#' @rdname internal

raise_yield_error <-
  function(model = "model",
           what = "property",
           try_this = NULL) {
    # assemble error message
    if (!is.null(try_this)) {
      try_mess <- paste0("\nTry ", try_this, ".")
    } else {
      try_mess = ""
    }
    mess <-
      paste0("Failed to yield the ", what, " from the ", model, ".", try_mess)
    rlang::abort(mess)
  }

#' @rdname internal
#' @keywords internal

get_response_from_formula <- function(x) {
  if (!rlang::is_formula(x)) {
    rlang::abort("x is not a formula.")
  }
  x_terms <- stats::terms(x)
  x_pos_response <- attr(x_terms, "response")
  if (x_pos_response == 0) {
    rlang::warn("x does not have a response")
    return(NULL)
  } else {
    return(rlang::as_string(as.list(attr(
      x_terms, "variables"
    ))[[x_pos_response + 1]]))
  }
}

#' @rdname internal
#' @keywords internal

abort_not_chantrics <- function(x) {
  if (!("chantrics" %in% class(x))) {
    rlang::abort("x is not a chantrics object", class = "chantrics_not_chantrics_object")
  }
}

#' @rdname internal
#' @keywords internal

get_variable_str_from_chantrics <- function(x) {
  abort_not_chantrics(x)
  return(paste(names(attr(x, "free_pars")), collapse = ", "))
}

#' @rdname internal
#' @keywords internal

get_formula_str_from_chantrics <- function(x) {
  abort_not_chantrics(x)
  return(paste(deparse(attr(x, "formula")), collapse = "\n"))
}
