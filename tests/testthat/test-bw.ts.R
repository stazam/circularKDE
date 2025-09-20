test_that("bw.ts returns a numeric value for valid input", {
  set.seed(60)
  x <- rvonmises(50, circular(pi), 1)
  result <- bw.ts(x)
  expect_equal(result, 2.301015)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("bw.ts returns a numeric value for valid input", {
  set.seed(123)
  x <- rvonmises(50, circular(pi / 2), 1)
  result <- bw.ts(x)
  expect_equal(result, 2.942911)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("bw.ts throws error on empty input", {
  expect_error(bw.ts(numeric(0)), "`x` must be a non-empty object.")
})

test_that("bw.ts throws error if input is not numeric", {
  expect_error(bw.ts(c("a", "b")), "must be a numeric vector")
})

test_that("bw.ts throws error if x contains only NAs", {
  x <- c(NA, NA)
  expect_error(bw.ts(x), "`x` contains all missing values.")
})

test_that("bw.ts removes NA values and returns result", {
  x <- circular(c(0, pi / 2, NA, pi))
  result <- bw.ts(x)
  expect_type(result, "double")
  expect_cli_warning(bw.ts(x),
                     1,
                     "! `x` contains missing values, which will be removed.")
})
