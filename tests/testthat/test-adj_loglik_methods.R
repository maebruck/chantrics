context("Tests concerning methods of adj_logLik()")

# ===== unit tests for adj_loglik() =====

test_that("adj_loglik adds the family name to the name of the returned object", {
  expect_equal(attr(fm_pois_adj, "name"), "poisson_glm_lm")
})

test_that("adj_loglik aborts if an invalid model type is given", {
  expect_error(adj_loglik("i'm a string"), class = "chantrics_invalid_model")
})

test_that("adj_loglik returns correct class", {
  expect_equal(class(fm_pois_adj)[1:2], c("chantrics", "chandwich"))
})

# definitions of fm_pois_adj objects in testthat.R

# ===== unit tests for anova.chantrics =====

test_that("anova.chantrics warns user if unnamed parameters are dropped", {
  expect_warning(anova(fm_pois_adj, fm_pois_small_adj, "foo", type = "vertical"),
    class = "chantrics_unnamed_params_dropped"
  )
})

test_that("anova.chantrics aborts if two models have the same # of params", {
  expect_error(anova(fm_pois_small_adj, adj_loglik(stats::update(
    fm_pois,
    formula = . ~ . - x
  ))), class = "chantrics_equal_num_params")
})

test_that("anova.chantrics aborts if attr(model, 'name') do not match", {
  expect_error(anova(fm_pois_adj, fm_negbin_adj), class = "chantrics_model_does_not_match")
})

test_that("anova.chantrics warns if the response changes", {
  fm_pois_small_adj_change_resp <- fm_pois_small_adj
  attr(fm_pois_small_adj_change_resp, "formula") <- z ~ x
  print(attr(fm_pois_adj, "formula"))
  print(attr(fm_pois_small_adj_change_resp, "formula"))
  expect_warning(anova(fm_pois_adj, fm_pois_small_adj_change_resp), class = "chantrics_response_param_discrepancy")
})

test_that("anova.chantrics aborts if the parameters are not a subset", {
  expect_error(anova(fm_pois_adj, fm_pois_cube_only_adj), class = "chantrics_params_not_subset")
})

test_that("anova.chantrics' sequential anova function returns the correct number of rows", {
  expect_equal(nrow(anova(fm_pois_adj)), 3)
})

# ==== nobs.chantrics() ====

test_that("nobs returns correct value", {
  expect_equal(stats::nobs(fm_pois_adj), stats::nobs(fm_pois))
})

# ==== update.chantrics() ====

test_that("update returns correct reduced model", {
  expect_equal(attr(update(fm_pois_adj, formula = . ~ . - I(x^2)), "adjSE"), attr(fm_pois_small_adj, "adjSE"))
})

test_that("update aborts when passed evaluate = FALSE", {
  expect_error(update(fm_pois_adj, evaluate = FALSE), class = "chantrics_update_evaluate_false")
})

# ==== alrtest() ====

test_that("alrtest executes properly", {
  expect_equal(nrow(alrtest(fm_pois_adj, fm_pois_small_adj, fm_pois_smallest_adj)), 3)
  expect_equal(nrow(alrtest(fm_pois_adj, 2, 1)), 3)
  expect_equal(nrow(alrtest(fm_pois_adj, "x", "I(x^2)")), 3)
  expect_equal(nrow(alrtest(fm_pois_adj, "x", 2)), 3)
  expect_equal(nrow(alrtest(fm_pois_adj, as.formula(". ~ . - I(x^2)"), as.formula(". ~ . - x"))), 3)
})

test_that("alrtest performs sequential anova if only passed one object", {
  expect_gt(nrow(alrtest(fm_pois_adj)), 1)
})

test_that("alrtest aborts if first object is not chantrics", {
  expect_error(alrtest("i'm character!"), class = "chantrics_not_chantrics_object")
})

test_that("alrtest aborts if second object is not supported", {
  expect_error(alrtest(fm_pois_adj, fm_pois), class = "chantrics_alrtest_unexpected_second_obj")
})

test_that("alrtest aborts if third objects are not of the right type", {
  expect_error(alrtest(fm_pois_adj, fm_pois_small_adj, "i'm character!"), class = "chantrics_not_chantrics_object")
  expect_error(alrtest(fm_pois_adj, 2, fm_pois), class = "chantrics_alrtest_numchar_match_failed")
  # expect_error(alrtest(fm_pois_adj, 2, 3.1415), class = "chantrics_alrtest_failed_int_coercion")
  expect_error(alrtest(fm_pois_adj, as.formula(". ~ . - x"), "dlr"), class = "chantrics_alrtest_not_formula")
})

test_that("alrtest aborts if the specified variable can't be found", {
  expect_error(alrtest(fm_pois_adj, 1456), class = "chantrics_alrtest_num_too_high")
  expect_error(alrtest(fm_pois_adj, "northernline"), class = "chantrics_alrtest_string_failed_matching")
})

# ==== df.residuals.chantrics() ====

test_that("df.residuals works as expected", {
  expect_equal(df.residual(fm_pois_adj), 247)
})

# ==== logLik_vec.chantrics() ====

test_that("logLik_vec.chantrics() components sum correctly", {
  expect_equal(sum(logLik_vec(fm_pois_adj)), as.numeric(logLik(fm_pois_adj)))
})

# ==== fitted.chantrics() ====

test_that("fitted.chantrics() generates values close to original ones", {
  expect_equal(stats::fitted(fm_pois_adj), stats::fitted(fm_pois))
})

# ==== residuals.chantrics() ====

test_that("residuals.chantrics() generates values close to original ones", {
  types <- c("response", "working", "pearson")
  for (type in types) {
    expect_equal(stats::residuals(fm_pois_adj, type = type), stats::residuals(fm_pois, type = type))
  }
})

# ==== predict.chantrics() ====

test_that("predict.chantrics() generates values close to original ones", {
  expect_equal(stats::predict(fm_pois_adj, type = "link"), stats::predict(fm_pois, type = "link"))
  expect_equal(stats::predict(fm_pois_adj, type = "response"), stats::predict(fm_pois, type = "response"))
  expect_equal(stats::predict(fm_negbin_theta_adj, type = "link"), stats::predict(fm_negbin_theta, type = "link"))
  expect_equal(stats::predict(fm_negbin_theta_adj, type = "response"), stats::predict(fm_negbin_theta, type = "response"))
})
