context("Tests logLik_vec() and its methods")

# logLik_vec() method tests go in the model's test file, not here.

# test logLik.logLik_vec() by comparing the result of the evaluation with a
# ("known") object's logLik (see glm)

# general functionality tests of the methods that plug onto logLik_vec()
llv_fm_pois <- logLik_vec(fm_pois)

test_that("logLik.logLik_vec() correctly sums all values together", {
  expect_equal(as.numeric(logLik(llv_fm_pois)), sum(llv_fm_pois))
})

test_that("nobs.logLik_vec() correctly returns the number of obserations", {
  expect_equal(nobs(llv_fm_pois), nobs(fm_pois))
})
