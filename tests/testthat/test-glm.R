context("Tests the aLogLik method and submethods associated with glm models.")

# poisson

# fit the misspecified poisson model from Introducing chandwich
# fitting happens in testthat.R
# also all definitions of fm_pois, fm_pois_adj, ...
# reference
pois_glm_loglik <- function(pars, y, x) {
  log_mu <- pars[1] + pars[2] * x + pars[3] * x^2
  return(dpois(y, lambda = exp(log_mu), log = TRUE))
}
reference_pois_logLik <- pois_glm_loglik(fm_pois$coefficients, y_nbinom, x_nbinom)
chantrics_pois_logLik <- logLik_vec(fm_pois, fm_pois$coefficients)

# ==== generate logistic regression data ====
# https://uvastatlab.github.io/2019/05/04/simulating-a-logistic-regression-model/
set.seed(1)
a_logit <- rbinom(250, 1, 0.5)
b_logit <- runif(250, 10, 80)
eta_logit <- 2 + 0.3 * a_logit - 0.02 * b_logit
# checked implementation with C source https://github.com/wch/r-source/blob/5a156a0865362bb8381dcd69ac335f5174a4f60c/src/library/stats/src/family.c#L73
probs_logit <- exp(eta_logit) / (1 + exp(eta_logit))
y_logit <- rbinom(n = 250, 1, probs_logit)
df_logit <- data.frame(y = y_logit, a = a_logit, b = b_logit)
logit_glm_loglik <- function(pars, df_logit) {
  eta <- pars[1] + pars[2] * df_logit$a + pars[3] * df_logit$b
  p <- exp(eta) / (1 + exp(eta))
  return(dbinom(df_logit$y, 1, p, log = TRUE))
}
bm_logit <- glm(y ~ a + b, data = df_logit, family = binomial(link = "logit"))
bm_logit_adj <- adj_loglik(bm_logit)
reference_logit_logLik <- logit_glm_loglik(bm_logit$coefficients, df_logit)
chantrics_logit_logLik <- logLik_vec(bm_logit, bm_logit$coefficients)


# ==== Generate Probit regression data ====
set.seed(1)
a_probit <- rbinom(250, 1, 0.5)
b_probit <- runif(250, 10, 80)
eta_probit <- 2 + 0.3 * a_probit - 0.02 * b_probit
# link is cdf of standard normal.
probs_probit <- pnorm(eta_probit, mean = 0, sd = 1)
y_probit <- rbinom(n = 250, 1, probs_probit)
df_probit <- data.frame(y = y_probit, a = a_probit, b = b_probit)
probit_glm_loglik <- function(pars, df_probit) {
  eta <- pars[1] + pars[2] * df_probit$a + pars[3] * df_probit$b
  p <- pnorm(eta, mean = 0, sd = 1)
  return(dbinom(df_probit$y, 1, p, log = TRUE))
}
bm_probit <- glm(y ~ a + b, data = df_probit, family = binomial(link = "probit"))
bm_probit_adj <- adj_loglik(bm_probit)
reference_probit_logLik <- probit_glm_loglik(bm_probit$coefficients, df_probit)
chantrics_probit_logLik <- logLik_vec(bm_probit, bm_probit$coefficients)
summary(bm_probit)
summary(bm_probit_adj)




test_that("logLik_vec.glm() returns correct loglik-vector if passed the correct object", {
  expect_equal(c(reference_pois_logLik), unname(c(chantrics_pois_logLik)))
  expect_equal(c(reference_logit_logLik), unname(c(chantrics_logit_logLik)))
  expect_equal(c(reference_probit_logLik), unname(c(chantrics_probit_logLik)))
})

test_that("logLik(logLik_vec.glm()) sums the log-likelihood correctly", {
  expect_equal(logLik(fm_pois), logLik(chantrics_pois_logLik))
  expect_equal(logLik(bm_logit), logLik(chantrics_logit_logLik))
  expect_equal(logLik(bm_probit), logLik(chantrics_probit_logLik))
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
  expect_error(model_generics_caller(bm_logit_adj))
})

## === ANOVA ===

test_that("Has the ANOVA function changed its output?", {
  expect_equivalent(round(anova(fm_pois_adj)[["I(x^2)", "ALRTS"]], 5), 1.82017)
  expect_equivalent(round(anova(fm_pois_adj, fm_pois_small_adj)[["2", "ALRTS"]], 5), 1.82017)
})
