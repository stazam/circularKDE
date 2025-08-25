test_that("bw.fo returns a numeric value for valid input", {
  set.seed(60)
  x <- rvonmises(50, circular(pi / 2), 1, control.circular = list(units = "radians"))
  result <- bw.fo(x)
  expect_equal(result, 0.228282328)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("bw.fo returns a numeric value for valid input", {
  set.seed(123)
  x <- rvonmises(50, circular(pi / 2), 1)
  result <- bw.fo(x)
  expect_equal(result, 0.598002222)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("bw.fo throws error on empty input", {
  expect_error(bw.fo(numeric(0)), "`x` must be a non-empty object.")
})

test_that("bw.fo throws error if input is not numeric", {
  expect_error(bw.fo(c("a", "b")), "must be a numeric vector")
})

test_that("bw.fo throws error if x contains only NAs", {
  x <- c(NA, NA)
  expect_error(bw.fo(x), "`x` contains all missing values.")
})

test_that("bw.fo removes NA values and returns result", {
  x <- circular(c(0, pi / 2, NA, pi))
  result <- bw.fo(x)
  expect_type(result, "double")
  expect_cli_warning(bw.fo(x),
                     1,
                     "! `x` contains missing values, which will be removed")
})

