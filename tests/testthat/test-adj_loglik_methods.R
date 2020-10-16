context("Tests concerning methods of adj_logLik()")

# ===== unit tests for adj_loglik() =====

test_that("adj_loglik adds the family name to the name of the returned object",
          {
            expect_equal(attr(fm_pois_adj, "name"), "poisson_glm_lm")
          })

test_that("adj_loglik aborts if an invalid model type is given", {
  expect_error(adj_loglik("i'm a string"), class = "chantrics_invalid_model")
})

test_that("adj_loglik returns correct class", {
  expect_equal(class(fm_pois_adj)[1:2], c("chantrics", "chandwich"))
})

#definitions of fm_pois_adj objects in testthat.R

# ===== unit tests for anova.chantrics =====

test_that("anova.chantrics warns user if unnamed parameters are dropped", {
  expect_warning(anova(fm_pois_adj, fm_pois_small_adj, "foo", type = "vertical"),
                 class = "chantrics_unnamed_params_dropped")
})

test_that("anova.chantrics aborts when no model is supplied", {
  expect_error(anova.chantrics("foo", "bar", type = "vertical"), class = "chantrics_not_enough_models")
})

test_that("anova.chantrics aborts if two models have the same # of params", {
  expect_error(anova(fm_pois_small_adj, adj_loglik(stats::update(
    fm_pois, formula = . ~ . - x
  ))), class = "chantrics_equal_num_params")
})

test_that("anova.chantrics aborts if attr(model, 'name') do not match", {
  expect_error(anova(fm_pois_adj, fm_negbin_small_adj), class = "chantrics_model_does_not_match")
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

# ==== nobs.chantrics() ====

test_that("nobs returns correct value", {
  expect_equal(stats::nobs(fm_pois_adj), stats::nobs(fm_pois))
})

# ==== update.chantrics() ====

test_that("update returns correct reduced model", {
  expect_equal(attr(update(fm_pois_adj, formula = . ~ . - I(x ^ 2)), "adjSE"), attr(fm_pois_small_adj, "adjSE"))
})

test_that("update aborts when passed evaluate = FALSE", {
  expect_error(update(fm_pois_adj, evaluate = FALSE), class = "chantrics_update_evaluate_false")
})
