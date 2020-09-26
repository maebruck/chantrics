#' Evaluate loglikelihood contributions from specific observations
#'
#' Generic function for calculating loglikelihood contributions from
#' individual observations for a fitted model.
#'
#' @param object A fitted model object.
#' @param ... Further arguments.
#' @export

logLik_vec <- function(object, ...) {
  UseMethod("logLik_vec")
}

#' Log-likelihood adjustments for fitted models
#'
#' This function adjusts the log-likelihood of fitted model objects based on
#' \href{http://doi.org/10.1093/biomet/asm015}{Chandler and Bate (2007)}.
#'
#' @export

# if required, turn this into a method (see logLik_vec) and the below into the
# .default() method

adj_loglik <- function(x,
                       cluster = NULL,
                       use_vcov = TRUE,
                       ...) {
  #check if x is a supported model type
  #adjust x
  adjusted_x <-
    chant_obj(x, cluster = cluster, use_vcov = use_vcov, ...)
  class(adjusted_x) <- c("chantrics", "chandwich", class(x))
  try(attr(adjusted_x, "formula") <-
        stats::formula(x), silent = TRUE)
  return(adjusted_x)
}
#' ANOVA tables: compare nested models
#'
#' \code{anova} method for \code{chantrics} objects
#'
#' Compare two or more nested models that have been adjusted using the
#' \code{\link{adj_logLik}} method. It uses the adjusted likelihood ratio test
#' statistic (ALRTS), as described in Section 3.5 of
#' \href{http://doi.org/10.1093/biomet/asm015}{Chandler and Bate (2007)}.
#'
#' Note that the objects specified do not have to be sorted in a specific way,
#' they will be sorted by the function by the number of parameters, as returned
#' by \code{attr(model1, "p_current")}.
#'
#' @param model1 Object of class \code{chantrics}, as returned by
#'   \code{\link{adj_logLik}}.
#' @param model2 Object of class \code{chantrics}, as returned by
#'   \code{\link{adj_logLik}}.
#' @param ... Further objects of class \code{chantrics}, as returned by
#'   \code{\link{adj_logLik}}, and/or parameters that will be passed to
#'   \code{\link[chandwich]{anova.chandwich}}, and then further to
#'   \code{\link[chandwich]{compare_models}}. The type of adjustment, out of
#'   \code{"vertical"}, \code{"cholesky"}, \code{"spectral"}, \code{"none"}, as
#'   specified in the parameter \code{type}, can also be specified here.
#'
#' @export
anova.chantrics <- function(model1, model2, ...) {
  dotargs <- list(...)
  potential_model_objects <-
    c(model1, model2, subset(dotargs, names(dotargs) == ""), use.names = FALSE)
  #save dotargs which are named
  named_dotargs <- subset(dotargs, names(dotargs) != "")
  #check whether the unnamed objects are actually chantrics objects
  pmo_is_chantrics <-
    vapply(potential_model_objects, function(x) {
      is(x, "chantrics")
    }, logical(1))
  #warn user that we drop supplied unnamed arguments that are not chantrics
  #objects
  if (!all(pmo_is_chantrics)) {
    rlang::warn(
      paste0(
        "One or more of the unnamed objects supplied ",
        "are not 'chantrics' objects.\nThey have been dropped."
      ),
      class = "chantrics_unnamed_params_dropped"
    )
  }
  #drop non-chantrics objects
  model_objects <-
    subset(potential_model_objects, pmo_is_chantrics)

  #check if there are at least two chantrics objects after dropping
  if (length(model_objects) < 2) {
    rlang::abort(
      paste0(
        "Less than two 'chantrics' objects have been supplied, ",
        "but we require at least two.\nPlease supply more models."
      ),
      class = "chantrics_not_enough_models"
    )
  }
  #sort models by # of parameters
  n_params <-
    vapply(model_objects, function(x) {
      attr(x, "p_current")
    }, numeric(1))
  #abort if # of params is identical for two models
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
  #order models
  model_objects <- model_objects[model_order]
  n_models <- length(model_objects)
  #compare_models requires single pairs
  largest_m <- model_objects[[1]]
  result_df <-
    data.frame(
      formula = c(NA),
      variable_string = c(get_variable_str_from_chantrics(largest_m)),
      resid_df = c(NA),
      resid_alrts = c(NA),
      df = c(NA),
      alrts = c(NA),
      p_value = c(NA)
    )
  try(result_df[["formula"]] <- c(get_formula_str_from_chantrics(largest_m)))
  for (i in 1:(n_models - 1)) {
    larger_m <- model_objects[[i]]
    smaller_m <- model_objects[[i + 1]]
    #check that the smaller model is indeed nested
    #check if they are the same model type
    if (attr(larger_m, "name") != attr(smaller_m, "name")) {
      rlang::warn(
        paste0(
          "The objects do not seem to stem from the same model.\n",
          "attr(model, 'name') do not match."
        ),
        class = "chantrics_model_does_not_match"
      )
    }
    #check if the response differs. This is only possible if original models
    #have a formula(x) function, otherwise continue
    if ("formula" %in% names(attributes(smaller_m))) {
      smaller_formula <- attr(smaller_m, "formula")
      #save for result_df
      result_df_nr_formula <- c(get_formula_str_from_chantrics(smaller_m))
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
      #if formula not available, write NA to result_df
      result_df_nr_formula <- c(NA)
    }

    #get the parameter names from the two models
    larger_free_params <- attr(larger_m, "free_pars")
    larger_free_params_names <- names(larger_free_params)
    smaller_free_params <- attr(smaller_m, "free_pars")
    smaller_free_params_names <- names(smaller_free_params)
    #check if the params of the smaller model are a subset of the larger.
    if (!all(smaller_free_params_names %in% larger_free_params_names)) {
      rlang::abort("The models are not nested: The parameters are not subsets.",
                   class = "chantrics_params_not_subset")
    }
    #find the fixed parameters
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
    print(result)
    #append to results data.frame
    result_df_nr <- data.frame(
      formula = result_df_nr_formula,
      variable_string = c(get_variable_str_from_chantrics(smaller_m)),
      resid_df = attr(smaller_m, "nobs"),
      resid_alrts = c(NA),
      df = c(NA),
      alrts = c(NA),
      p_value = c(NA)
    )
  }
  #get names from objects
  #if all functions have formulae available, then get formulae


  heading <- c("Analysis of Adjusted Deviance Table\n")
  structure(table, heading = )
}
