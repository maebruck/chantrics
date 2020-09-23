library(testthat)
library(chantrics)

#function to check if all generic functions for fitted models are accessible

#how to run this test:
# test_that("Are generics accessible for adjusted glm models?", {
#   #test that there is no error
#   expect_error(model_generics_caller(), regexp = NA)
# })

model_generics_caller <- function(object) {
  print(object)
  summary(object)
  #coef() and coefficients() call the same S3 methods
  coef(object)
  #residuals() and resid() call the same S3 methods
  #residuals(object)
  #fitted(object)
  #anova() tested in test-anova.R
  #predict(object)
  #plot() requires a single free covariate
  if (attr(object, "p_current") == 1) {
    plot(object, type = 1:4)
  }
  confint(object)
  #deviance(object)
  vcov(object)
  logLik(object)
  AIC(object)
}

#define test object using glm and poisson
set.seed(123)
x <- rnorm(250)
y <- rnbinom(250, mu = exp(1 + x), size = 1)
fm_pois <- glm(y ~ x + I(x ^ 2), family = poisson)
fm_pois_adj <- adj_loglik(fm_pois)
fm_pois_small <- update(fm_pois, formula = . ~ . - I(x ^ 2))
fm_pois_small_adj <- adj_loglik(fm_pois_small)
fm_pois_smallest <- update(fm_pois, formula = . ~ 1)
fm_pois_smallest_adj <- adj_loglik(fm_pois_smallest)
fm_pois_cube_only <- update(fm_pois, formula = . ~ I(x ^ 3))
fm_pois_cube_only_adj <- adj_loglik(fm_pois_cube_only)
#just a dummy until implemented
fm_negbin_small_adj <- fm_pois_small_adj
attr(fm_negbin_small_adj, "name") <- "negbin_glm_lm"

test_check("chantrics")
