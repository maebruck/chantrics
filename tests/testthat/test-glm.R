context("Tests the aLogLik method and submethods associated with glm models.")

# fit the misspecified poisson model from Introducing chandwich
# fitting happens in testthat.R
# also all definitions of fm_pois, fm_pois_adj, ...
# reference
pois_glm_loglik <- function(pars, y, x) {
  log_mu <- pars[1] + pars[2] * x + pars[3] * x^2
  return(dpois(y, lambda = exp(log_mu), log = TRUE))
}
reference_pois_logLik <- pois_glm_loglik(fm_pois$coefficients, y, x)
chantrics_pois_logLik <- logLik_vec(fm_pois, fm_pois$coefficients)
test_that("logLik_vec.glm() returns correct loglik-vector if passed the correct object", {
  expect_equal(c(reference_pois_logLik), unname(c(chantrics_pois_logLik)))
})

test_that("logLik(logLik_vec.glm()) sums the log-likelihood correctly", {
  expect_equal(logLik(fm_pois), logLik(chantrics_pois_logLik))
  # add calculations for other families here
})

test_that("adj_logLik can handle use_vcov = F (set very high error tolerance.)", {
  expect_equal(summary(fm_pois_adj), summary(adj_logLik(fm_pois, use_vcov = F)), tolerance = 1e-3)
})

# !! add unit tests for other glm link functions for logLik_vec() here !!

# adjust fm_pois

test_that("Are generics accessible for adjusted glm models?", {
  # test that there is no error
  expect_error(model_generics_caller(fm_pois_adj), regexp = NA)
})

## === ANOVA ===

test_that("Has the ANOVA function changed its output?", {
  expect_equivalent(round(anova(fm_pois_adj)[["I(x^2)", "ALRTS"]], 5), 1.82017)
  expect_equivalent(round(anova(fm_pois_adj, fm_pois_small_adj)[["2", "ALRTS"]], 5), 1.82017)
})
