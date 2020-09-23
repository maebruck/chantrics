context("Tests concerning methods of adj_logLik()")


test_that("adj_loglik adds the family name to the name of the returned object",
          {
            expect_equal(attr(fm_pois_adj, "name"), "poisson_glm_lm")
          })


# !! add unit tests for adj_loglik here !!


#definitions of fm_pois_adj objects in testthat.R

# ===== unit tests for anova.chantrics =====

test_that("anova.chantrics warns user if unnamed parameters are dropped", {
  expect_warning(anova(fm_pois_adj, fm_pois_small_adj, "foo", e = "e"), class = "chantrics_unnamed_params_dropped")
})

test_that("anova.chantrics aborts when less than 2 models are supplied", {
  expect_error(anova(fm_pois_adj, "foo", "bar", e = "e"), class = "chantrics_not_enough_models")
})

test_that("anova.chantrics aborts if two models have the same # of params", {
  expect_error(anova(fm_pois_small_adj, adj_loglik(update(
    fm_pois, formula = . ~ . - x
  ))), class = "chantrics_equal_num_params")
})

test_that("anova.chantrics warns if attr(model, 'name') do not match", {
  expect_warning(anova(fm_pois_adj, fm_negbin_small_adj), class = "chantrics_model_does_not_match")
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
