#' chandwich internal functions
#'
#' Internal functions: These are not to be exported.
#' @name internal
#' @keywords internal
NULL

#' @keywords internal
#' @rdname internal

chant_obj <- function(x, cluster = NULL, ...){
  #add tests if object valid
  #does x have a logLik_vec method?
  if (!any(paste0("logLik_vec.", class(x)) %in% utils::methods("logLik_vec"))){
    rlang::abort("x does not have a logLik_vec method")
  }
  #create function for log-likelihood of x
  logLik_f <- function(pars, fitted_object, ...){
    return(c(logLik_vec(fitted_object, pars = pars)))
  }
  #get mle estimate from x
  mle = stats::coef(x)
  #adjust object using chandwich
  adjusted_obj <- chandwich::adjust_loglik(loglik = logLik_f, cluster = cluster, fitted_object = x, p = length(mle), par_names = names(mle), name = paste(class(x), collapse = "_", mle = mle))
  class(adjusted_obj) <- c("chantrics", "chandwich")
  return(adjusted_obj)
}

#' @keywords internal
#' @rdname internal
is.error <- function(x) inherits(x, "try-error")

#' @keywords internal
#' @rdname internal

raise_yield_error <- function(model = "model", what = "property", try_this = NULL){
  # assemble error message
  if (!is.null(try_this)){
    try_mess <- paste0("\nTry ", try_this, ".")
  } else {
    try_mess = ""
  }
  mess <- paste0("Failed to yield the ", what, " from the ", model, ".", try_mess)
  rlang::abort(mess)
}
