context("Tests the aLogLik method and submethods associated with glm models.")


# stuff about dispersion parameters
# https://stats.stackexchange.com/questions/33432/dispersion-parameter-in-glm-output

test_that("logLik_vec.glm() returns correct loglik-vector if passed the correct object", {
  expect_equal(c(reference_pois_logLik), unname(c(chantrics_pois_logLik)))
  expect_equal(c(reference_logit_logLik), unname(c(chantrics_logit_logLik)))
  expect_equal(c(reference_probit_logLik), unname(c(chantrics_probit_logLik)))
  expect_equal(c(reference_gauss_logLik), unname(c(chantrics_gauss_logLik)), tolerance = 1e-3)
  expect_equal(c(reference_negbin_logLik), unname(c(chantrics_negbin_logLik)))
  expect_equal(c(reference_negbin_theta_logLik), unname(c(chantrics_negbin_theta_logLik)))
})

test_that("logLik(logLik_vec.glm()) sums the loglikelihood correctly", {
  expect_equal(logLik(fm_pois), logLik(chantrics_pois_logLik))
  expect_equal(logLik(bm_logit), logLik(chantrics_logit_logLik))
  expect_equal(logLik(bm_probit), logLik(chantrics_probit_logLik))
  # expect_equal(logLik(glm_gauss), logLik(chantrics_gauss_logLik), tolerance = 1e-2)
  expect_equal(logLik(fm_negbin), logLik(chantrics_negbin_logLik))
  expect_equal(logLik(fm_negbin_theta), logLik(chantrics_negbin_theta_logLik))
  # add calculations for other families here
})

test_that("adj_logLik can handle use_vcov = F (set very high error tolerance.)", {
  expect_equal(summary(fm_pois_adj), summary(adj_loglik(fm_pois, use_vcov = F)), tolerance = 1e-3)
  expect_equal(summary(bm_logit_adj), summary(adj_loglik(bm_logit, use_vcov = F)), tolerance = 1e-3)
})

test_that("Are generics accessible for adjusted glm models?", {
  # test that there is no error
  expect_error(model_generics_caller(fm_pois_adj), regexp = NA)
  expect_error(model_generics_caller(bm_logit_adj), regexp = NA)
  expect_error(model_generics_caller(bm_probit_adj), regexp = NA)
  expect_error(model_generics_caller(glm_gauss_adj), regexp = NA)
  expect_error(model_generics_caller(fm_negbin_adj), regexp = NA)
  expect_error(model_generics_caller(fm_negbin_theta_adj, run.anova = TRUE), regexp = NA)
})

test_that("logLik_vec.glm() aborts/warns if pars does not conform to the design matrix", {
  pars <- c(1, 2, 3)
  names(pars) <- c("(Intercept)", "x", "I(x^2)")
  expect_error(logLik_vec(fm_pois, pars = pars[1:2]), class = "chantrics_pars_wrong_length")
  paes <- pars
  names(paes)[1] <- "(Interbebd)"
  expect_warning(logLik_vec(fm_pois, pars = paes), class = "chantrics_parnames_do_not_match")
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

## === dispersion_stat() ===

# test_that("Does dispersion.stat() calculate the correct dispersion parameter?", {
#   realtheta <- summary(fm_negbin_theta)$theta
#   testtheta <- chantrics:::dispersion.stat(y_nbinom, fitted(fm_negbin_theta), fm_negbin_theta)
#   expect_equal(testtheta, realtheta)
# })
