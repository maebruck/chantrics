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
#'   * [glm]
#'
#' @return An object of class `"chantrics"` inheriting from class `"chandwich"`.
#'   See [chandwich::adjust_loglik()]. The remaining elements of the returned
#'   class are `class(x)`.
#'
#'   `chantrics` objects have `AIC`, `anova`, `coef`, `confint`, `logLik`,
#'   `nobs`, `plot`, `print`, `summary` and `vcov` methods.
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
  supported_models <- c("glm")
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
  try(name_pieces <- c(x$family$family, name_pieces), silent = TRUE)
  # get mle estimate from x
  mle <- stats::coef(x)
  # estimate Hessian if use_vcov is TRUE
  if (use_vcov) {
    # only possible if vcov exists
    if (any(paste0("vcov.", class(x)) %in% utils::methods("vcov"))) {
      H <- -solve(stats::vcov(x))
    } else {
      # otherwise use integrated methods
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
  class(adjusted_x) <- c("chantrics", "chandwich", class(x))
  try(attr(adjusted_x, "formula") <-
    stats::formula(x), silent = TRUE)
  attr(adjusted_x, "unadj_object") <- x
  args_list <- rlang::dots_list(...)
  args_list[["cluster"]] <- cluster
  args_list[["use_vcov"]] <- use_vcov
  attr(adjusted_x, "chantrics_args") <- args_list
  return(adjusted_x)
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
    c(object, subset(dotargs, names(dotargs) == ""), use.names = FALSE)
  # save dotargs which are named
  named_dotargs <- subset(dotargs, names(dotargs) != "")
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
    # unadjusted_object <- attr(adjusted_object, "unadj_object")
    # get list of variables by which anova should split
    formula_full <-
      try(attr(prev_adjusted_object, "formula"), silent = TRUE)
    if (is.error(formula_full)) {
      rlang::abort(
        paste0(
          "Formula not available in object ",
          "via attr(model1, 'formula')\n",
          "Handling of this case not supported."
        )
      )
    }
    variable_vec <-
      rev(attr(stats::terms(formula_full), "term.labels"))
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
        "Less than two 'chantrics' objects have been supplied, ",
        "but we require at least two.\nPlease supply more models."
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
  result_df.resid_df[[1]] <- get_resid_df_from_chantrics(largest_m)
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
    # code below doesn't work
    # fixed_pars <- larger_free_params[fixed_pars]
    # fixed_at <- 0
    # print(fixed_pars)
    # result <- do.call(chandwich::compare_models, c(list(larger = larger_m, smaller = smaller_m, fixed_pars = fixed_pars, fixed_at = fixed_at), named_dotargs))
    result <-
      do.call(chandwich::compare_models, c(list(
        larger = larger_m, smaller = smaller_m
      ), named_dotargs))
    # append to results data.frame
    result_df.formula[[i + 1]] <- result_df_nr_formula
    result_df.resid_df[[i + 1]] <-
      get_resid_df_from_chantrics(smaller_m)
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
  orig_obj <- attr(object, "unadj_object")
  chantrics_args[["x"]] <-
    stats::update(orig_obj, ..., evaluate = TRUE)
  return(do.call(adj_loglik, chantrics_args))
}
