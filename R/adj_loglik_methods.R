#' Log-likelihood adjustments for fitted models
#'
#' This function adjusts the log-likelihood of fitted model objects based on
#' [Chandler and Bate (2007)](http://doi.org/10.1093/biomet/asm015). It is a
#' generic function for different types of models, which are listed in
#' **Supported models**. This section also contains links to function-specific
#' help pages.
#'
#' @param x A fitted model object that is supported, see **Supported models**
#'
#' @param cluster A vector or factor indicating from which cluster the
#'   respective log-likelihood contributions originate. Must have the same
#'   length as the vector returned by [logLik_vec()]. If `cluster` is not
#'   supplied or `NULL` then it is assumed that each observation forms its own
#'   cluster.
#'
#' @param use_vcov A logical scalar. By default, the [vcov()] method for `x` is
#'   used to estimate the Hessian of the independence loglikelihood, if the
#'   function exists. Otherwise, or if `use_vcov = FALSE`, `H` is estimated
#'   using [stats::optimHess()] inside [chandwich::adjust_loglik()].
#'
#' @param ... Further arguments to be passed to [sandwich::meatCL()] if
#'   `cluster` is defined, if `cluster = NULL`, they are passed into
#'   [sandwich::meat()].
#'
#' @details If `use_vcov = TRUE`, the current default, the function will test
#'   whether a `vcov` S3 method exists for `x`, and will take the
#'   variance-covariance matrix from there. Otherwise, or if `use_vcov = FALSE`
#'   the variance-covariance matrix of the MLE is estimated inside
#'   [chandwich::adjust_loglik()] using [stats::optimHess()].
#'
#' @section Supported models:
#' * [glm]
#' * [glm.nb]
#'
#' @return An object of class `"chantrics"` inheriting from class `"chandwich"`.
#'   See [chandwich::adjust_loglik()]. The remaining elements of the returned
#'   class are `class(x)`.
#'
#'   `chantrics` objects have `AIC`, `anova`, `coef`, `confint`, `df.residual`,
#'   `logLik`, `logLik_vec`, `nobs`, `plot`, `print`, `summary` and `vcov`
#'   methods.
#'
#' @section Examples: See the model-specific pages in the *supported models*
#'   section.
#'
#' @references R. Chandler and S. Bate, Inference for clustered data using the
#'   independence loglikelihood, Biometrika, 94 (2007), pp. 167–183.
#'   <http://doi.org/10.1093/biomet/asm015>.
#'
#' @seealso [lax::alogLik()] supports adjustment for user-supplied objects.
#'
#' @export

adj_loglik <- function(x,
                       cluster = NULL,
                       use_vcov = TRUE,
                       ...) {
  # if required, turn this into a method (see logLik_vec) and the below into the
  # .default() method
  supported_models <- c("glm", "negbin")
  # check if x is a supported model type
  if (!(class(x)[1] %in% supported_models)) {
    rlang::abort(
      paste0(
        class(x)[1],
        " is not a supported model type.\n",
        "Please refer to ?adj_loglik() for a list of valid model types."
      ),
      class = "chantrics_invalid_model"
    )
  }
  # adjust x
  if (!any(paste0("logLik_vec.", class(x)) %in% utils::methods("logLik_vec"))) {
    rlang::abort("x does not have a logLik_vec method")
  }
  # create function for log-likelihood of x
  logLik_f <- function(pars, fitted_object, ...) {
    return(c(logLik_vec(fitted_object, pars = pars)))
  }
  name_pieces <- c(class(x))
  # add glm family to name
  if (class(x)[1] == "glm") {
    try(
      {
        name_pieces <- c(x$family$family, name_pieces)
      },
      silent = TRUE
    )
  }
  # get mle estimate from x
  mle <- stats::coef(x)

  # estimate Hessian if use_vcov is TRUE
  if (use_vcov) {
    # only possible if vcov exists
    if (any(paste0("vcov.", class(x)) %in% utils::methods("vcov"))) {
      H <- -solve(stats::vcov(x))
    } else {
      # otherwise use integrated methods - emergency fallback.
      H <- NULL
    }
  } else {
    H <- NULL
  }
  if (is.null(cluster)) {
    V <-
      sandwich::meat(x, fitted_object = x, loglik_fn = logLik_f, ...) * stats::nobs(x)
  } else {
    V <-
      sandwich::meatCL(
        x,
        cluster = cluster,
        fitted_object = x,
        loglik_fn = logLik_f,
        ...
      ) * stats::nobs(x)
  }
  # adjust object using chandwich
  adjusted_x <-
    chandwich::adjust_loglik(
      loglik = logLik_f,
      fitted_object = x,
      p = length(mle),
      par_names = names(mle),
      name = paste(name_pieces, collapse = "_"),
      mle = mle,
      H = H,
      V = V
    )
  # post-estimation
  if (class(x)[1] == "glm") {
    if (x$family$family == "gaussian") {
      # estimate dispersion
      response_vec <- get_response_from_model(x)
      eta_vec <- get_design_matrix_from_model(x) %*% attr(adjusted_x, "res_MLE")
      mu_vec <- x$family$linkinv(eta_vec)
      attr(adjusted_x, "dispersion") <- dispersion.gauss(response_vec, mu_vec, stats::df.residual(x) - 1)
    } else if (substr(x$family$family, 1, 18) == "Negative Binomial(") {
      # estimate dispersion
      response_vec <- get_response_from_model(x)
      eta_vec <- get_design_matrix_from_model(x) %*% attr(adjusted_x, "res_MLE")
      mu_vec <- x$family$linkinv(eta_vec)
      attr(adjusted_x, "dispersion") <- dispersion.stat(response_vec, mu_vec, x)
    }
  } else if (class(x)[1] == "negbin") {
    # estimate dispersion
    response_vec <- get_response_from_model(x)
    eta_vec <- get_design_matrix_from_model(x) %*% attr(adjusted_x, "res_MLE")
    mu_vec <- x$family$linkinv(eta_vec)
    attr(adjusted_x, "theta") <- MASS::theta.ml(y = response_vec, mu = mu_vec)
    # hijack x to pass adjusted theta to original model object to circumvent reestimation for confint and anova
    attr(attr(adjusted_x, "loglik_args")[["fitted_object"]], "theta_chantrics") <- attr(adjusted_x, "theta")
    # attr(adjusted_x, "dispersion") <- dispersion.stat(response_vec, mu_vec, x)
  }
  class(adjusted_x) <- c("chantrics", "chandwich", class(x))
  try(attr(adjusted_x, "formula") <-
    stats::formula(x), silent = TRUE)
  args_list <- rlang::dots_list(...)
  args_list[["cluster"]] <- cluster
  args_list[["use_vcov"]] <- use_vcov
  attr(adjusted_x, "chantrics_args") <- args_list
  # check if unadjusted model has been passed through, if not, add it
  try(if (is.null(attr(adjusted_x, "loglik_args")[["fitted_object"]])) {
    attr(adjusted_x, "loglik_args")[["fitted_object"]] <- x
  })
  return(adjusted_x)
}

#' @export

summary.chantrics <- function(object, ...) {
  # https://stackoverflow.com/a/8316856/
  ans <- NextMethod()
  if (!is.null(attr(object, "dispersion"))) {
    attr(ans, "dispersion") <- attr(object, "dispersion")
  }
  if (!is.null(attr(object, "theta"))) {
    attr(ans, "theta") <- attr(object, "theta")
  }
  class(ans) <- c("summary.chantrics", class(ans))
  return(ans)
}

#' @export

print.summary.chantrics <- function(x, digits = max(3L, getOption("digits") - 3L),
                                    # symbolic.cor = x$symbolic.cor,
                                    signif.stars = getOption("show.signif.stars"), ...) {
  stats::printCoefmat(x, digits = digits)
  if (!is.null(attr(x, "theta"))) {
    cat("\n(Theta: ", format(attr(x, "theta"), digits = digits), ", SE: ", format(attr(attr(x, "theta"), "SE"), digits = digits), ")\n", sep = "")
  }
  if (!is.null(attr(x, "dispersion"))) {
    cat("\n(Dispersion parameter taken to be ", format(attr(x, "dispersion"), digits = digits), ")\n", sep = "")
  }
  invisible(x)
}


#' ANOVA tables: compare nested models
#'
#' `anova` method for `chantrics` objects
#'
#' Create an analysis of adjusted deviance table for one object (sequential), or
#' two or more nested models that have been adjusted using the
#' [adj_loglik()] method. It uses the adjusted likelihood ratio test
#' statistic (ALRTS), as described in Section 3.5 of
#' [Chandler and Bate (2007)](http://doi.org/10.1093/biomet/asm015).
#'
#' @param object Object of class `chantrics`, as returned by
#'   [adj_loglik()].
#' @param ... Further objects of class `chantrics`, as returned by
#'   [adj_loglik()], and/or parameters that will be passed to
#'   [chandwich::compare_models()]. The type of adjustment, out of
#'   `"vertical"`, `"cholesky"`, `"spectral"`, `"none"`, as
#'   specified in the parameter `type`, can also be specified here.
#'
#' @details Each line represents the model as given above the table, with each
#'   line (except for the first line) showing the residual degrees of freedom of
#'   that model, and the change in degrees of freedom, the ALRTS and the
#'   associated p-value in comparison to the model in the line above.
#'
#'   When a single model is specified, the function returns a sequential
#'   analysis of deviance table, where, iteratively, one term is being removed
#'   from the right of the full formula. This process is continued until the
#'   "intercept only" model is left. The row names are the names of the dropped
#'   term in comparison to the model in the line above.
#'
#'   If more than one model is specified, the function sorts the models by their
#'   number of variables as returned by [adj_loglik()] in
#'   `attr(x, "p_current")`.

#'
#' Details of the ALRT can be found in [chandwich::compare_models()] and in
#' [Chandler and Bate (2007)](http://doi.org/10.1093/biomet/asm015).
#'
#' @return An object of class `"anova"` inheriting from class `"data.frame"`.
#'   The columns are as follows: \item{Resid.df}{The residual number of degrees
#'   of freedom in the model.} \item{df}{The increase in residual degrees of
#'   freedom with respect to the model in the row above.} \item{ALRTS}{The
#'   adjusted likelihood ratio statistic.} \item{Pr(>ALRTS)}{The p-value of the
#'   test that the model above is a "significantly better" model as the one in
#'   the current row.}
#'
#' @references R. Chandler and S. Bate, Inference for clustered data using the
#'   independence loglikelihood, Biometrika, 94 (2007), pp.
#'   167–183. <http://doi.org/10.1093/biomet/asm015>.
#'
#' @seealso [chandwich::compare_models]: implementation of the comparison
#'   mechanism
#'
#' @seealso [chandwich::anova.chandwich]: `anova` method of the `chandwich`
#'   package, which also uses `compare.models()`
#'
#' @examples
#'
#' # from Introducing Chandwich.
#' set.seed(123)
#' x <- rnorm(250)
#' y <- rnbinom(250, mu = exp(1 + x), size = 1)
#' fm_pois <- glm(y ~ x + I(x^2), family = poisson)
#' fm_pois_adj <- adj_loglik(fm_pois)
#' fm_pois_small <- update(fm_pois, formula = . ~ . - I(x^2))
#' fm_pois_small_adj <- adj_loglik(fm_pois_small)
#' fm_pois_smallest <- update(fm_pois, formula = . ~ 1)
#' fm_pois_smallest_adj <- adj_loglik(fm_pois_smallest)
#'
#' anova(fm_pois_adj, fm_pois_small_adj, fm_pois_smallest_adj)
#' # use different types of adjustment with type, default is "vertical"
#' anova(fm_pois_adj, fm_pois_small_adj, fm_pois_smallest_adj, type = "cholesky")
#'
#' # sequential anova
#' anova(fm_pois_adj)
#' @export
#'

## S3 method for class 'chantrics'
anova.chantrics <- function(object, ...) {
  dotargs <- rlang::dots_list(...)
  potential_model_objects <-
    c(object, get_unnamed(dotargs), use.names = FALSE)
  # save dotargs which are named
  named_dotargs <- get_named(dotargs)
  # check whether the unnamed objects are actually chantrics objects
  pmo_is_chantrics <-
    vapply(potential_model_objects, function(x) {
      is(x, "chantrics")
    }, logical(1))
  # warn user that we drop supplied unnamed arguments that are not chantrics
  # objects
  if (!all(pmo_is_chantrics)) {
    rlang::warn(
      paste0(
        "One or more of the unnamed objects supplied ",
        "are not 'chantrics' objects.\nThey have been dropped."
      ),
      class = "chantrics_unnamed_params_dropped"
    )
  }
  # drop non-chantrics objects
  model_objects <-
    subset(potential_model_objects, pmo_is_chantrics)

  # check if there is at least one chantrics objects after dropping
  if (length(model_objects) == 1) {
    # ==== Sequential ANOVA ====
    # create sequential model objects
    prev_adjusted_object <- model_objects[[1]]
    # unadjusted_object <- attr(adjusted_object, "loglik_args")[["fitted_object"]]
    # get list of variables by which anova should split
    variable_vec <-
      rev(attr(stats::terms(prev_adjusted_object), "term.labels"))
    # initialise progress bar
    pb <- progress::progress_bar$new(total = length(variable_vec))
    for (rm_this_var in variable_vec) {
      pb$tick()
      prev_adjusted_object <-
        stats::update(prev_adjusted_object, formula = stats::as.formula(paste0(". ~ . - ", rm_this_var)))
      model_objects <- c(model_objects, prev_adjusted_object)
    }
  } else if (length(model_objects) < 1) {
    rlang::abort(
      paste0(
        "No 'chantrics' objects have been supplied, ",
        "but we require at least one.\nPlease supply more models."
      ),
      class = "chantrics_not_enough_models"
    )
  }
  # sort models by # of parameters
  n_params <-
    vapply(model_objects, function(x) {
      attr(x, "p_current")
    }, numeric(1))
  # abort if # of params is identical for two models
  if (anyDuplicated(n_params)) {
    rlang::abort(
      paste0(
        "More than one supplied chantrics model ",
        "has the same number of parameters.\n",
        "Please ensure that all models are nested with decreasing ",
        "parameter count."
      ),
      class = "chantrics_equal_num_params"
    )
  }
  model_order <- order(n_params, decreasing = TRUE)
  # order models
  model_objects <- model_objects[model_order]
  n_models <- length(model_objects)
  # compare_models requires single pairs
  largest_m <- model_objects[[1]]

  # create vector for each result_df entry, and join at the end, is more
  # efficient, see
  # https://stackoverflow.com/questions/20689650/how-to-append-rows-to-an-r-data-frame
  result_df.formula <- character(n_models)
  result_df.formula[[1]] <-
    get_variable_str_from_chantrics(largest_m)
  try(result_df.formula[[1]] <-
    get_formula_str_from_chantrics(largest_m),
  silent = TRUE
  )
  result_df.resid_df <- integer(n_models)
  result_df.resid_df[[1]] <- df.residual.chantrics(largest_m)
  result_df.df <- integer(n_models)
  result_df.df[[1]] <- NA
  result_df.alrts <- numeric(n_models)
  result_df.alrts[[1]] <- NA
  result_df.p_value <- numeric(n_models)
  result_df.p_value[[1]] <- NA

  for (i in 1:(n_models - 1)) {
    larger_m <- model_objects[[i]]
    smaller_m <- model_objects[[i + 1]]
    # check that the smaller model is indeed nested
    # check if they are the same model type
    if (attr(larger_m, "name") != attr(smaller_m, "name")) {
      # chandwich::compare_models does the same check, and aborts, -> abort
      # here too.
      rlang::abort(
        paste0(
          "The objects do not seem to stem from the same model.\n",
          "attr(model, 'name') do not match."
        ),
        class = "chantrics_model_does_not_match"
      )
    }
    # check if the response differs. This is only possible if original models
    # have a formula(x) function, otherwise continue
    if ("formula" %in% names(attributes(smaller_m))) {
      smaller_formula <- attr(smaller_m, "formula")
      # save for result_df
      result_df_nr_formula <-
        get_formula_str_from_chantrics(smaller_m)
      if ("formula" %in% names(attributes(larger_m))) {
        larger_response <-
          get_response_from_formula(attr(larger_m, "formula"))
        smaller_response <-
          get_response_from_formula(attr(smaller_m, "formula"))
        if (larger_response != smaller_response) {
          rlang::warn(
            paste0(
              "The response parameter in attr(model, 'formula')",
              "do not seem to match."
            ),
            class = "chantrics_response_param_discrepancy"
          )
        }
      }
    } else {
      # if formula not available, write variable list to result_df
      # emergency fallback
      result_df_nr_formula <-
        get_variable_str_from_chantrics(smaller_m)
    }

    # get the parameter names from the two models
    larger_free_params <- attr(larger_m, "free_pars")
    larger_free_params_names <- names(larger_free_params)
    smaller_free_params <- attr(smaller_m, "free_pars")
    smaller_free_params_names <- names(smaller_free_params)
    # check if the params of the smaller model are a subset of the larger.
    if (!all(smaller_free_params_names %in% larger_free_params_names)) {
      rlang::abort("The models are not nested: The parameters are not subsets.",
        class = "chantrics_params_not_subset"
      )
    }
    # find the fixed parameters
    index_fixed_pars <-
      which(!(larger_free_params_names %in% smaller_free_params_names))
    fixed_pars <- larger_free_params[index_fixed_pars]
    fixed_at <- rep(0, length(fixed_pars))
    names(fixed_at) <- names(fixed_pars)
    attr(larger_m, "fixed_pars") <- NULL
    attr(larger_m, "fixed_at") <- NULL
    attr(smaller_m, "fixed_pars") <- fixed_pars
    attr(smaller_m, "fixed_at") <- fixed_at
    result <-
      do.call(chandwich::compare_models, c(list(
        larger = larger_m, smaller = smaller_m
      ), named_dotargs))
    # append to results data.frame
    result_df.formula[[i + 1]] <- result_df_nr_formula
    result_df.resid_df[[i + 1]] <- df.residual.chantrics(smaller_m)
    result_df.df[[i + 1]] <- result[["df"]]
    result_df.alrts[[i + 1]] <- result[["alrts"]]
    result_df.p_value[[i + 1]] <- result[["p_value"]]
  }
  # create data.frame from vectors
  result_df <- data.frame(
    resid_df = result_df.resid_df,
    df = result_df.df,
    alrts = result_df.alrts,
    p_value = result_df.p_value
  )
  dimnames(result_df)[[2]] <-
    c("Resid.df", "df", "ALRTS", "Pr(>ALRTS)")
  try(dimnames(result_df)[[1]] <-
    c("full model", variable_vec),
  silent = TRUE
  )
  title <- "Analysis of Adjusted Deviance Table\n"
  topnote <-
    paste0("Model ",
      format(1:n_models),
      ": ",
      result_df.formula,
      collapse = "\n"
    )
  structure(
    result_df,
    heading = c(title, topnote),
    class = c("anova", "data.frame")
  )
}


#' @importFrom stats nobs
#' @export
#'

nobs.chantrics <- function(object, ...) {
  if (!missing(...)) {
    rlang::warn("extra arguments discarded")
  }
  return(attr(object, "nobs"))
}

#' Update, re-fit and re-adjust a Model Call
#'
#' `update.chantrics()` will update a model that has been adjusted by
#' [adj_loglik()]. It passes all arguments to the standard [stats::update()]
#' function.
#'
#' The function will not update any arguments passed to the `adj_loglik()`
#' function, re-run `adj_loglik()` with the changed arguments.
#'
#' @param object A `"chantrics"` returned by [adj_loglik()].
#' @param ... Additional arguments to the call, passed to [stats::update()] to
#'   update the original model specification.
#'
#' @details Passing `evaluate = FALSE` is not supported, if this is required,
#'   run [stats::update()] on the unadjusted object.
#'
#' @return The fitted, adjusted `"chantrics"` object.
#'
#' @seealso [stats::update()]
#' @seealso [stats::update.formula()]
#'
#' @examples
#' # from Introducing Chandwich.
#' set.seed(123)
#' x <- rnorm(250)
#' y <- rnbinom(250, mu = exp(1 + x), size = 1)
#' fm_pois <- glm(y ~ x + I(x^2), family = poisson)
#' fm_pois_adj <- adj_loglik(fm_pois)
#' fm_pois_small_adj <- update(fm_pois_adj, formula = . ~ . - I(x^2))
#' summary(fm_pois_small_adj)
#' fm_pois_smallest_adj <- update(fm_pois_adj, formula = . ~ 1)
#' summary(fm_pois_smallest_adj)
#' @importFrom stats update
#' @export

update.chantrics <- function(object, ...) {
  # break if evaluate is false
  dotargs <- rlang::dots_list(...)
  if ("evaluate" %in% names(dotargs)) {
    if (!dotargs[["evaluate"]]) {
      rlang::abort("'evaluate = FALSE' is not supported.", class = "chantrics_update_evaluate_false")
    }
  }
  chantrics_args <- get_additional_args_from_chantrics_call(object)
  orig_obj <- get_unadj_object(object)
  chantrics_args[["x"]] <-
    stats::update(orig_obj, ..., evaluate = TRUE)
  return(do.call(adj_loglik, chantrics_args))
}

#' @importFrom stats terms
#' @export

terms.chantrics <- function(x, ...) {
  # passes terms object through from unadjusted object
  return(terms(attr(x, "loglik_args")[["fitted_object"]]))
}

#' Adjusted Likelihood Ratio Test of Nested Models
#'
#' `alrtest` is a helper function to simulate the functions
#' [lmtest::waldtest()]/[lmtest::lrtest()] for adjusted `chantrics` objects. The
#' method can be employed to compare nested models (see details).
#'
#' @param object a `chantrics` object as returned from [adj_loglik()].
#' @param ... further object specifications and/or parameters that will be
#'   passed to [chandwich::compare_models()]. The type of adjustment, out of
#'   `"vertical"`, `"cholesky"`, `"spectral"`, `"none"`, as specified in the
#'   parameter `type`, can also be specified here.
#'
#' @details This function is a helper function that creates an interface to
#'   [anova.chantrics()] that is similar to [lmtest::waldtest()] and
#'   [lmtest::lrtest()].
#'
#'   The standard method is to compare the fitted model object `object` with the
#'   models in `...`. Instead of passing the fitted models into `...`, other
#'   specifications are possible. Note that the types of specifications can't be
#'   mixed, except between numerics/characters. The type of the second object
#'   supplied determines the algorithm used.
#'   * **`"chantrics"` objects**: When
#'   supplying two or more `"chantrics"` objects, they will be sorted as in
#'   [anova.chantrics()]. Then, the ALRTS will be computed consecutively between
#'   the two neighbouring models. Note that all models must be nested. For
#'   details refer to [anova.chantrics()].
#'   * **`"numeric"`**: If the second
#'   object is `"numeric"` or `"character"`, then `"numeric"` objects
#'   corresponding element in `attr(terms(object1), "term.labels")` will be
#'   turned into their corresponding `"character"` element and will be handled
#'   as in `"character"` below.
#'   * **`"character"`**: If the second object is
#'   `"numeric"` or `"character"`, then the `"character"` objects are
#'   consecutively included in an update formula like `update(object1, . ~ . -
#'   object2)`
#'   * **`"formula"`**: If the second object is a `"formula"`, then
#'   the model will be computed as `update(object1, object2)`.
#'
#'   Then, the adjusted likelihood ratio test statistic (ALRTS), as described in
#'   Section 3.5 of [Chandler and Bate
#'   (2007)](http://doi.org/10.1093/biomet/asm015), is computed by
#'   [anova.chantrics()].
#'
#'   If a single unnamed object is passed in `...`, sequential ANOVA is perfomed
#'   on `object`.
#'
#' @return An object of class `"anova"` inheriting from class `"data.frame"`.
#'   The columns are as follows: \item{Resid.df}{The residual number of degrees
#'   of freedom in the model.} \item{df}{The increase in residual degrees of
#'   freedom with respect to the model in the row above.} \item{ALRTS}{The
#'   adjusted likelihood ratio statistic.} \item{Pr(>ALRTS)}{The p-value of the
#'   test that the model above is a "significantly better" model as the one in
#'   the current row.}
#'
#' @seealso [anova.chantrics()] for the implementation of the computations of
#'   the test statistics.
#' @seealso [lmtest::waldtest()] and [lmtest::lrtest()] for syntax.
#'
#' @export

alrtest <- function(object, ...) {
  # object has to be chantrics.
  abort_not_chantrics(object)
  dotargs <- rlang::dots_list(...)
  # store named arguments to be passed down into the anova() call
  named_args <- get_named(dotargs)
  unnamed_args <- unname(get_unnamed(dotargs))
  # create object_list that is object + unnamed_args
  if (length(unnamed_args) == 0) {
    # if only one chantrics object is passed in, run sequential anova.
    return(do.call("anova", c(object, named_args)))
  }
  object_list <- c(object, unnamed_args)

  # functions that handle the conversion of the unnamed args into chantrics objects.
  # note that alrtest will not check if what's done to the models actually makes sense
  # this is left to anova.chantrics.
  chantrics_handler <- function(object_list) {
    # check if all objects are of type chantrics, since anova only warns if model objects don't fit
    checks_logi <- vapply(object_list, function(x) {
      class(x)[[1]] != "chantrics"
    }, logical(1))
    if (any(checks_logi)) {
      index_false <- as.character(which(checks_logi))
      rlang::abort(
        paste0(
          "Objects in position ",
          index_false,
          "are not chantrics objects.\n",
          "Since the second object is chantrics, all other unnamed objects must be chantrics."
        ),
        class = "chantrics_not_chantrics_object"
      )
    }
    return(object_list)
  }
  numchar_handler <- function(object_list) {
    chantrics_obj <- object_list[[1]]
    # flatten any structures, e.g. lists or vectors into one big list
    flat_object_list <- purrr::flatten(object_list[-1])
    # get a list of covariate names from the first object
    var_names <- attr(terms(chantrics_obj), "term.labels")
    # function to pass into vapply
    match_var_names <- function(x) {
      # check that all elements in flat_object_list are numericals or characters
      if (rlang::is_scalar_double(x)) {
        # try coercing into integer - emergency fallback.
        x_int <- try(as.integer(x), silent = TRUE)
        if (is.error(x_int)) {
          rlang::abort("Could not coerce ",
            as.character(x),
            "into an integer.",
            class = "chantrics_alrtest_failed_int_coercion"
          )
        }
        # try to match the numerical elements with elements in var_names
        if (x_int > length(var_names)) {
          rlang::abort(
            paste0(
              as.character(x_int),
              " exceeds ",
              length(var_names),
              ", the number of variables in\n",
              "attr(terms(object1), 'term.labels')"
            ),
            class = "chantrics_alrtest_num_too_high"
          )
        }
        return(var_names[[x_int]])
      } else if (rlang::is_string(x)) {
        # try to match the string against the list of variable names, abort if
        # not matched
        if (x %in% var_names) {
          return(x)
        } else {
          rlang::abort(
            paste0(
              x,
              " could not be matched against any element in\n",
              "attr(terms(object1), 'term.labels')"
            ),
            class = "chantrics_alrtest_string_failed_matching"
          )
        }
      } else {
        rlang::abort(
          paste0(
            "The object ",
            as.character(x),
            " of class ",
            class(x),
            "\n",
            "is neither coercible into an integer nor a matchable character string."
          ),
          class = "chantrics_alrtest_numchar_match_failed"
        )
      }
    }
    remove_these_vars <-
      vapply(flat_object_list, match_var_names, character(1))
    ready_objects <- list(chantrics_obj)
    # then, pass the objects one by one into update functions, consecutively
    # decreasing the list of variables leave error handling and matching errors
    # to update().
    for (remove_var in remove_these_vars) {
      ready_objects <-
        c(
          ready_objects,
          stats::update(ready_objects[[length(ready_objects)]], formula = stats::as.formula(paste0(
            ". ~ . - ", remove_var
          )))
        )
    }
    return(ready_objects)
  }
  formula_handler <- function(object_list) {
    ready_objects <- list(object_list[[1]])
    # pass the formula objects consecutively into update. Leave all error handling to formula.
    for (curr_formula in object_list[-1]) {
      # check if formula
      if (!rlang::is_formula(curr_formula)) {
        rlang::abort(paste0(as.character(curr_formula), " is not a formula."), class = "chantrics_alrtest_not_formula")
      } else {
        ready_objects <-
          c(
            ready_objects,
            stats::update(ready_objects[[length(ready_objects)]], formula = curr_formula)
          )
      }
    }
    return(ready_objects)
  }
  # write bank of if statements that check if the second element in object_list
  # is a chantrics, a numerical/character, or a formula object then lead the
  # object_list into the handlers
  checkobj <- object_list[[2]]
  if (class(checkobj)[[1]] == "chantrics") {
    ready_objects <- chantrics_handler(object_list)
  } else if (any(is.character(checkobj), is.numeric(checkobj))) {
    ready_objects <- numchar_handler(object_list)
  } else if (rlang::is_formula(checkobj)) {
    ready_objects <- formula_handler(object_list)
  } else {
    # if it is none of the types, abort
    rlang::abort(paste0("The second object is of unexpected class", class(checkobj)),
      class = "chantrics_alrtest_unexpected_second_obj"
    )
  }
  # pass ready_objects together with named_args into anova.chantrics
  anova_result <-
    do.call("anova", c(unname(ready_objects), named_args))
  attr(anova_result, "heading")[[1]] <-
    "Adjusted likelihood ratio test\n"
  return(anova_result)
}

#' @importFrom stats df.residual
#' @export

df.residual.chantrics <- function(object, ...) {
  abort_not_chantrics(object)
  if (!missing(...)) {
    rlang::warn("extra arguments discarded")
  }
  return(stats::nobs(object) - attr(object, "p_current"))
}

#' @export

logLik_vec.chantrics <- function(object, ...) {
  if (!missing(...)) {
    rlang::warn("extra arguments discarded")
  }
  return(attr(object, "loglikVecMLE"))
}

#' @importFrom stats confint
#' @export

confint.chantrics <- function(object, ...) {
  if (!is.null(attr(object, "theta"))) {
    # save theta in external env to discourage re-calculation
    # https://stackoverflow.com/a/10904331/
    assign("theta", attr(object, "theta"), envir = bypasses.env)
  }
  res <- NextMethod()
  # http://adv-r.had.co.nz/Environments.html#env-basics
  if (exists("theta", envir = bypasses.env, inherits = FALSE)) {
    rm("theta", envir = bypasses.env)
  }
  return(res)
}

#' @importFrom stats fitted
#' @export

fitted.chantrics <- function(object, ...) {
  if (!missing(...)) {
    rlang::warn("extra arguments discarded")
  }
  abort_not_chantrics(object)
  modelname <- unlist(strsplit(attr(object, "name"), "_"))
  if ("glm" %in% modelname) {
    fittedvals <- fittedhelper.glm(object, modelname)
  } else {
    rlang::abort(paste0("'",attr(object, "name"),"' is currently not supported."))
  }
  return(fittedvals)
}

#' @importFrom stats residuals
#' @export

residuals.chantrics <- function(object, type = c("response", "working", "pearson"), ...) {
  # https://www.datascienceblog.net/post/machine-learning/interpreting_generalized_linear_models/
  if (!missing(...)) {
    rlang::warn("extra arguments discarded")
  }
  abort_not_chantrics(object)
  type <- match.arg(type)
  modelname <- unlist(strsplit(attr(object, "name"), "_"))
  if ("glm" %in% modelname) {
    response <- get_response_from_model(object)
    fitted_responselev <- fitted(object)
    resids_responselev <- response - fitted_responselev
  } else {
    rlang::abort(paste0("'",attr(object, "name"),"' is currently not supported."))
  }
  if (type == "response") {
    # standard residuals on response level
    result <- resids_responselev
  } else if (type == "working") {
    # response residuals normalised by fitted values
    result <- resids_responselev / fitted_responselev
  } else if (type == "pearson") {
    # response residuals normalised by sqrt of the estimate
    result <- resids_responselev / sqrt(fitted_responselev)
  } else {
    rlang::abort(paste0("Residuals of type '", type, "' are not implemented."))
  }
  return(result)
}


#' Predict Method for chantrics fits
#'
#' Obtains predictions from chantrics objects.
#'
#' @param object Object of class `chantrics`, as returned by `adj_loglik()`
#' @param newdata optionally, a data frame in which to look for variables with which to predict. If omitted, the fitted linear predictors are used.
#' @param type the type of prediction required. The default `"response"` is on the scale of the response variables. The alternative `"link"` is on the scale of the linear predictors, if applicable, otherwise, an error is returned.
#'
#' @details If `newdata` is omitted, the predictions are based on the data used for the fit.
#' Any instances of `NA` will return `NA`.
#'
#' @return A vector of predictions.
#'
#' @importFrom stats predict
#'

predict.chantrics <- function(object, newdata = NULL, type = c("response", "link")) {
  type <- match.arg(type)
  # get the data frame with all observations
  if (missing(newdata)) {
    newdata <- get_design_matrix_from_model(object)
  }
  # check that the required parameters are available
  if (inherits(object, "glm")) {
    # get list of required coefficients
    coef_names <- names(stats::coef(object))
    # match these with newdata and create new matrix with only those columns
    design_matrix <- newdata[coef_names]
    eta_vec <- object %*% attr(object, "res_MLE")
    mu_vec <- object$family$linkinv(eta_vec)
  }
}
