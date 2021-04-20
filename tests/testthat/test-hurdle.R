context("Tests the aLogLik method and submethods associated with glm models.")

test_that("adj_loglik() aborts if the model matrices are not accessible", {
  rd_hurdle_no_mm <- pscl::hurdle(trips ~ . | quality + income, data = RecreationDemand, dist = "negbin", x = FALSE)
  expect_error(adj_loglik(rd_hurdle_no_mm), class = "chantrics_missing_model_matrices")
})

# test_that("logLik_vec.glm() returns correct loglik-vector if passed the correct object", {
#   #expect_equal(c(reference_negbin_theta_logLik), unname(c(chantrics_negbin_theta_logLik)))
# })

# test_that("logLik(logLik_vec.glm()) sums the loglikelihood correctly", {
#   #expect_equal(logLik(fm_pois), logLik(chantrics_pois_logLik))
# })

# test_that("adj_logLik can handle use_vcov = F (set very high error tolerance.)", {
#   #expect_equal(summary(fm_pois_adj), summary(adj_loglik(fm_pois, use_vcov = F)), tolerance = 1e-3)
# })

test_that("Are generics accessible for adjusted hurdle models?", {
  # test that there is no error
  expect_error(model_generics_caller(rd_hurdle_nb_adj, run.anova = FALSE, run.confint = FALSE, run.residuals = FALSE), regexp = NA)
  expect_error(model_generics_caller(rd_hurdle_logit_poi_adj, run.anova = FALSE, run.confint = FALSE, run.residuals = FALSE), regexp = NA)
  expect_error(model_generics_caller(rd_hurdle_geom_geom_adj, run.anova = FALSE, run.confint = FALSE, run.residuals = FALSE), regexp = NA)
})
