test_that("bw.jf returns a numeric value for valid input", {
  set.seed(60)
  x <- rvonmises(50, circular(pi / 2), 1, control.circular = list(units = "radians"))
  result <- bw.jf(x)
  expect_equal(result, 1.38316346)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("bw.jf returns a numeric value for valid input", {
  set.seed(123)
  x <- rvonmises(50, circular(pi / 2), 1)
  result <- bw.jf(x)
  expect_equal(result, 1.78506578)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("bw.jf throws error on empty input", {
  expect_error(bw.jf(numeric(0)), "`x` must be a non-empty object.")
})

test_that("bw.jf throws error if input is not numeric", {
  expect_error(bw.jf(c("a", "b")), "must be a numeric vector")
})

test_that("bw.jf throws error if x contains only NAs", {
  x <- c(NA, NA)
  expect_error(bw.jf(x), "`x` contains all missing values.")
})

test_that("bw.jf removes NA values and returns result", {
  x <- circular(c(0, pi / 2, NA, pi))
  result <- bw.jf(x)
  expect_type(result, "double")
  expect_cli_warning(bw.jf(x),
                     1,
                     "! `x` contains missing values, which will be removed.")
})

