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

test_check("chantrics")
