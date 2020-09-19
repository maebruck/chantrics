context("Tests the aLogLik method and submethods associated with glm models.")

test_that("logLik_vec.glm() returns correct loglik-vector if passed the correct object",{
  #from Introducing Chandwich
  set.seed(123)
  x <- rnorm(250)
  y <- rnbinom(250, mu = exp(1 + x), size = 1)
  fm_pois <- glm(y ~ x + I(x^2), family = poisson)
  #reference
  pois_glm_loglik <- function(pars, y, x) {
    log_mu <- pars[1] + pars[2] * x + pars[3] * x ^ 2
    return(dpois(y, lambda = exp(log_mu), log = TRUE))
  }
  expect_equal(c(pois_glm_loglik(fm_pois$coefficients, y, x)), unname(c(logLik_vec(fm_pois, fm_pois$coefficients))))
})

