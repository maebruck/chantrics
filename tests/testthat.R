library(testthat)
library(chantrics)

# function to check if all generic functions for fitted models are accessible

# how to run this test:
# test_that("Are generics accessible for adjusted glm models?", {
#   #test that there is no error
#   expect_error(model_generics_caller(), regexp = NA)
# })

model_generics_caller <- function(object, run.anova = TRUE, run.confint = TRUE, run.residuals = TRUE) {
  invisible(print(object))
  summary(object)
  # coef() and coefficients() call the same S3 methods
  coef(object)
  # residuals() and resid() call the same S3 methods
  if (run.residuals) {
    residuals(object)
    fitted(object)
  }
  # predict(object)
  # plot() requires a single free covariate
  if (attr(object, "p_current") == 1) {
    invisible(plot(object, type = 1:4))
  }
  if (run.confint) {
    confinttest <- confint(object)
    plot(confinttest)
  }
  # deviance(object)
  vcov(object)
  logLik(object)
  AIC(object)
  # sequential anova?
  if (run.anova) {
    anova(object)
  }
  df.residual(object)
  nobs(object)
}

# ===== glm test models =====

# ==== generate negbin data ====

# fit the misspecified poisson model from Introducing chandwich

set.seed(123)
x_nbinom <- rnorm(250)
y_nbinom <- rnbinom(250, mu = exp(1 + x_nbinom), size = 1)
df_nbinom <- data.frame(x = x_nbinom, y = y_nbinom)

# ==== poisson test objects ====

fm_pois <- glm(y ~ x + I(x^2), data = df_nbinom, family = poisson)
fm_pois_adj <- adj_loglik(fm_pois)
fm_pois_small <- update(fm_pois, formula = . ~ . - I(x^2))
fm_pois_small_adj <- adj_loglik(fm_pois_small)
fm_pois_smallest <- update(fm_pois, formula = . ~ 1)
fm_pois_smallest_adj <- adj_loglik(fm_pois_smallest)
fm_pois_cube_only <- update(fm_pois, formula = . ~ I(x^3))
fm_pois_cube_only_adj <- adj_loglik(fm_pois_cube_only)

pois_glm_loglik <- function(pars, y, x) {
  log_mu <- pars[1] + pars[2] * x + pars[3] * x^2
  return(dpois(y, lambda = exp(log_mu), log = TRUE))
}
reference_pois_logLik <- pois_glm_loglik(fm_pois$coefficients, y_nbinom, x_nbinom)
chantrics_pois_logLik <- logLik_vec(fm_pois, fm_pois$coefficients)

# ==== negbin test objects ====
if (!requireNamespace("MASS", quietly = TRUE)) {
  rlang::abort("requires MASS")
}
fm_negbin <- glm(y ~ x, data = df_nbinom, family = MASS::negative.binomial(theta = 1))
# summary(fm_negbin)
# summary(fm_pois_small)
fm_negbin_adj <- adj_loglik(fm_negbin)
# summary(fm_negbin_adj)
negbin_glm_loglik <- function(pars, df, theta = 1) {
  eta <- pars[1] + pars[2] * df$x
  return(dnbinom(df$y, size = theta, mu = exp(eta), log = TRUE))
}
reference_negbin_logLik <- negbin_glm_loglik(fm_negbin$coefficients, df_nbinom, theta = 1)
chantrics_negbin_logLik <- logLik_vec(fm_negbin, fm_negbin$coefficients)



# estimation of theta
fm_negbin_theta <- MASS::glm.nb(y ~ x, data = df_nbinom)
# summary(fm_negbin_theta)
fm_negbin_theta_adj <- adj_loglik(fm_negbin_theta)
# summary(fm_negbin_theta_adj)
reference_negbin_theta_logLik <- negbin_glm_loglik(fm_negbin_theta$coefficients, df_nbinom, theta = summary(fm_negbin_theta)$theta)
chantrics_negbin_theta_logLik <- chantrics:::logLik_vec(fm_negbin_theta, fm_negbin_theta$coefficients)


# ==== generate gaussian data ====
# following example in https://www.r-bloggers.com/2013/08/fitting-a-model-by-maximum-likelihood/

set.seed(1001)
sample_gaussian <- 250
a_gauss <- rbinom(sample_gaussian, 1, 0.5)
b_gauss <- runif(sample_gaussian, 10, 80)
y_gauss <- 2 + 0.3 * a_gauss + 0.95 * b_gauss + rnorm(sample_gaussian, mean = 0, sd = 2)
df_gauss <- data.frame(y = y_gauss, a = a_gauss, b = b_gauss)
glm_gauss <- glm(y ~ a + b, data = df_gauss, family = gaussian())
glm_gauss_adj <- adj_loglik(glm_gauss)
gauss_glm_loglik <- function(pars, disp, df) {
  eta <- df$y - pars[1] - pars[2] * df$a - pars[3] * df$b
  return(dnorm(eta, 0, sd = sqrt(disp), log = TRUE))
}
reference_gauss_logLik <- gauss_glm_loglik(glm_gauss$coefficients, summary(glm_gauss)$dispersion, df_gauss)
chantrics_gauss_logLik <- logLik_vec(glm_gauss, glm_gauss$coefficients)


# ==== generate logistic regression data ====
# https://uvastatlab.github.io/2019/05/04/simulating-a-logistic-regression-model/
set.seed(1)
sample_logit <- 250
a_logit <- rbinom(sample_logit, 1, 0.5)
b_logit <- runif(sample_logit, 10, 80)
eta_logit <- 2 + 0.3 * a_logit - 0.02 * b_logit
# checked implementation with C source https://github.com/wch/r-source/blob/5a156a0865362bb8381dcd69ac335f5174a4f60c/src/library/stats/src/family.c#L73
probs_logit <- exp(eta_logit) / (1 + exp(eta_logit))
y_logit <- rbinom(n = sample_logit, 1, probs_logit)
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
sample_probit <- 250
a_probit <- rbinom(sample_probit, 1, 0.5)
b_probit <- runif(sample_probit, 10, 80)
eta_probit <- 2 + 0.3 * a_probit - 0.02 * b_probit
# link is cdf of standard normal.
probs_probit <- pnorm(eta_probit, mean = 0, sd = 1)
y_probit <- rbinom(n = sample_probit, 1, probs_probit)
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

# ===== hurdle models =====

# ==== Estimate hurdle model ====

if (!requireNamespace("pscl", quietly = TRUE)) {
  rlang::abort("requires pscl")
}
data("RecreationDemand", package = "AER")
rd_hurdle_nb <- pscl::hurdle(trips ~ . | quality + income, data = RecreationDemand, dist = "negbin", x = TRUE)
rd_hurdle_nb_adj <- adj_loglik(rd_hurdle_nb)
rd_hurdle_nb_small <- pscl::hurdle(trips ~ . - income | quality + income, data = RecreationDemand, dist = "negbin", x = TRUE)
rd_hurdle_nb_small_adj <- adj_loglik(rd_hurdle_nb_small)

rd_hurdle_logit_poi <- pscl::hurdle(trips ~ . | quality + income, data = RecreationDemand, dist = "poisson", zero.dist = "binomial", link = "logit", x = TRUE)
rd_hurdle_logit_poi_adj <- adj_loglik(rd_hurdle_logit_poi)

rd_hurdle_geom_geom <- pscl::hurdle(trips ~ . | quality + income, data = RecreationDemand, dist = "geometric", zero.dist = "geometric", x = TRUE)
rd_hurdle_geom_geom_adj <- adj_loglik(rd_hurdle_geom_geom)
#


test_check("chantrics")

