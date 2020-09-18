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
