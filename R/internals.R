#' chandwich internal functions
#'
#' Internal functions: These are not to be exported.
#' @name internal
#' @keywords internal
NULL

#' @keywords internal
#' @rdname internal

chand_obj <- function(x, cluster = NULL, ...){
  #add tests if object valid

  #adjust object using chandwich
  adjusted_obj <- chandwich::adjust_loglik(loglik = )
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
  abort(mess)
}
