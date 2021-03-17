#' chantrics: Loglikelihood Adjustments for Econometric Models
#'
#' `chantrics` adjusts the loglikelihood of common econometric models for
#' clustered data based on the estimation process suggested in
#' [Chandler and Bate (2007)](http://doi.org/10.1093/biomet/asm015),
#' using the [chandwich](https://cran.r-project.org/package=chandwich) package,
#' and provides convenience functions for inference on the adjusted models.
#' `adj_loglik()` adjusts the parameter covariance matrix of the models to
#' incorporate clustered data, and can mitigate model misspecification by
#' wrapping `chandwich::adjust_loglik` for the supported models.
#'
#' The returned model of class `chantrics` can be plugged into standard model
#' evaluation and model comparison methods, for example `summary()`, `confint()`
#' and `anova()`, and a hypothesis test framework provided by `alrtest()`.
#'
#' See \code{vignette("chantrics-vignette", package = "chantrics")} for an
#' overview of the package.
#'
#' @references R. Chandler and S. Bate, Inference for clustered data using the
#'   independence loglikelihood, Biometrika, 94 (2007), pp. 167–183.
#'   <http://doi.org/10.1093/biomet/asm015>.
#'
#' @docType package
#' @name chantrics
NULL
