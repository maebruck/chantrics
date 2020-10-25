#' chandwich internal functions
#'
#' Internal functions: These are not to be exported.
#' @name internal
#' @keywords internal
NULL

#' @keywords internal
#' @rdname internal
is.error <- function(x) {
  inherits(x, "try-error")
}

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
      try_mess <- ""
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

#' @rdname internal
#' @keywords internal

get_resid_df_from_chantrics <- function(x) {
  abort_not_chantrics(x)
  return(stats::nobs(x) - attr(x, "p_current"))
}

#' @rdname internal
#' @keywords internal

get_additional_args_from_chantrics_call <- function(x) {
  abort_not_chantrics(x)
  # parse call and remove function name and model arg.
  # return(as.list(str2lang(attr(x, "chantrics_call")))[-c(1, 2)])
  return(attr(x, "chantrics_args"))
}

#' @rdname internal
#' @keywords internal

get_named <- function(x){
  return(subset(x, names(x) != ""))
}

#' @rdname internal
#' @keywords internal

get_unnamed <- function(x){
  return(subset(x, names(x) == ""))
}
