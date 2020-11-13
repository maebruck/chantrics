context("Tests the aLogLik method and submethods associated with glm models.")


#stuff about dispersion parameters
#https://stats.stackexchange.com/questions/33432/dispersion-parameter-in-glm-output

test_that("logLik_vec.glm() returns correct loglik-vector if passed the correct object", {
  expect_equal(c(reference_pois_logLik), unname(c(chantrics_pois_logLik)))
  expect_equal(c(reference_logit_logLik), unname(c(chantrics_logit_logLik)))
  expect_equal(c(reference_probit_logLik), unname(c(chantrics_probit_logLik)))
  expect_equal(c(reference_gauss_logLik), unname(c(chantrics_gauss_logLik)))
})

test_that("logLik(logLik_vec.glm()) sums the log-likelihood correctly", {
  expect_equal(logLik(fm_pois), logLik(chantrics_pois_logLik))
  expect_equal(logLik(bm_logit), logLik(chantrics_logit_logLik))
  expect_equal(logLik(bm_probit), logLik(chantrics_probit_logLik))
  expect_equal(logLik(glm_gauss), logLik(chantrics_gauss_logLik))
  # add calculations for other families here
})

test_that("adj_logLik can handle use_vcov = F (set very high error tolerance.)", {
  expect_equal(summary(fm_pois_adj), summary(adj_loglik(fm_pois, use_vcov = F)), tolerance = 1e-3)
  expect_equal(summary(bm_logit_adj), summary(adj_loglik(bm_logit, use_vcov = F)), tolerance = 1e-3)
})

# !! add unit tests for other glm link functions for logLik_vec() here !!

# adjust fm_pois

test_that("Are generics accessible for adjusted glm models?", {
  # test that there is no error
  expect_error(model_generics_caller(fm_pois_adj), regexp = NA)
  expect_error(model_generics_caller(bm_logit_adj), regexp = NA)
  expect_error(model_generics_caller(bm_probit_adj), regexp = NA)
  expect_error(model_generics_caller(glm_gauss_adj), regexp = NA)
})

## === ANOVA ===

test_that("Has the ANOVA function changed its output?", {
  expect_equivalent(round(anova(fm_pois_adj)[["I(x^2)", "ALRTS"]], 5), 1.82017)
  expect_equivalent(round(anova(fm_pois_adj, fm_pois_small_adj)[["2", "ALRTS"]], 5), 1.82017)
})

## === dispersion.gauss() ===

test_that("Does dispersion.gauss() calculate the correct dispersion parameter?", {
  realdisp <- summary(glm_gauss)$dispersion
  testdisp <- chantrics:::dispersion.gauss(y_gauss, fitted(glm_gauss), stats::nobs(glm_gauss) - 3)
  expect_equal(realdisp, testdisp)
})
