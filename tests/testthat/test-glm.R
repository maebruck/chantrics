context("Tests the aLogLik method and submethods associated with glm models.")

#fit the misspecified poisson model from Introducing chandwich
set.seed(123)
x <- rnorm(250)
y <- rnbinom(250, mu = exp(1 + x), size = 1)
fm_pois <- glm(y ~ x + I(x ^ 2), family = poisson)
#reference
pois_glm_loglik <- function(pars, y, x) {
  log_mu <- pars[1] + pars[2] * x + pars[3] * x ^ 2
  return(dpois(y, lambda = exp(log_mu), log = TRUE))
}
reference_pois_logLik <- pois_glm_loglik(fm_pois$coefficients, y, x)
chantrics_pois_logLik <- logLik_vec(fm_pois, fm_pois$coefficients)
test_that("logLik_vec.glm() returns correct loglik-vector if passed the correct object",
          {expect_equal(c(reference_pois_logLik), unname(c(chantrics_pois_logLik)))
          })

# !! add unit tests for other glm link functions for logLik_vec() here !!

#adjust fm_pois
fm_pois_adj <- adj_loglik(fm_pois)

# !! add unit tests for adj_loglik here !!

test_that("Are generics accessible for adjusted glm models?", {
  #test that there is no error
  expect_error(model_generics_caller(fm_pois_adj), regexp = NA)
})

