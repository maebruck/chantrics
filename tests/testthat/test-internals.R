context("Tests internal functions")

# is.error()

test_that("is.error() identifies a try-error object as true, and nothing else", {
  expect_true(is.error(try(log("a"), silent = TRUE)))
  expect_false(is.error("string"))
  expect_false(is.error(10))
})


# raise_yield_error()

test_that("raise_yield_error(a,b,c) correctly outputs an error", {
  expect_error(raise_yield_error())
})

# abort_not_chantrics()

test_that("abort_not_chantrics(x) aborts if x is not a chantrics object", {
  expect_error(abort_not_chantrics("foo"), class = "chantrics_not_chantrics_object")
  expect_error(abort_not_chantrics(fm_pois_adj), regexp = NA)
})

# get_variable_str_from_chantrics()

test_that("get_variable_str_from_chantrics aborts if x is not a chantrics object", {
  expect_error(get_variable_str_from_chantrics("foo"), class = "chantrics_not_chantrics_object")
})

test_that("get_variable_str_from_chantrics returns a character string with the variable names", {
  expect_equal(get_variable_str_from_chantrics(fm_pois_adj), "(Intercept), x, I(x^2)")
})

# get_formula_str_from_chantrics

test_that("get_formula_str_from_chantrics aborts if x is not a chantrics object", {
  expect_error(get_formula_str_from_chantrics("foo"), class = "chantrics_not_chantrics_object")
})

test_that("get_formula_str_from_chantrics returns a character string with the formula", {
  expect_equal(get_formula_str_from_chantrics(fm_pois_adj), "y ~ x + I(x^2)")
})
