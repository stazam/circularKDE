test_that("bw.lscvg returns a numeric value for valid input", {
  set.seed(60)
  x <- rvonmises(50, circular(pi / 2), 1, control.circular = list(units = "radians"))
  result <- bw.lscvg(x)
  expect_equal(result, 6.6326176)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("bw.lscvg returns a numeric value for valid input", {
  set.seed(123)
  x <- rvonmises(50, circular(pi / 2), 1)
  result <- bw.lscvg(x)
  expect_equal(result, 5.09811023)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("bw.lscvg throws error on empty input", {
  expect_error(bw.lscvg(numeric(0)), "`x` must be a non-empty object.")
})

test_that("bw.lscvg throws error if input is not numeric", {
  expect_error(bw.lscvg(c("a", "b")), "must be a numeric vector")
})

test_that("bw.lscvg throws error if x contains only NAs", {
  x <- c(NA, NA)
  expect_error(bw.lscvg(x), "`x` contains all missing values.")
})

test_that("bw.lscvg removes NA values and returns result", {
  x <- c(0, pi / 2, NA, pi)
  result <- bw.lscvg(x)
  expect_type(result, "double")
  expect_cli_warning(bw.lscvg(x),
                     1,
                     "! `x` contains missing values, which will be removed")
})

test_that("bw.lscvg handles non-numeric g", {
  x <- seq(0, 2 * pi, length.out = 5)
  expect_cli_warning(
    result <- bw.lscvg(x, g = "wrong"),
    1,
    "! Argument `g` must be numeric. Default value 4 for coefficient was used."
  )
  expect_type(result, "double")
})

test_that("bw.lscvg handles non-numeric lower", {
  x <- seq(0, 2 * pi, length.out = 5)
  expect_cli_warning(
    result <- bw.lscvg(x, lower = "zero"),
    1,
    "! Argument `lower` must be numeric. Default value 0 for lower boundary was used."
  )
  expect_type(result, "double")
})

test_that("bw.lscvg handles non-numeric upper", {
  x <- seq(0, 2 * pi, length.out = 5)
  expect_cli_warning(
    result <- bw.lscvg(x, upper = "sixty"),
    1,
    "! Argument `upper` must be numeric. Default value 60 for upper boundary was used."
  )
  expect_type(result, "double")
})

test_that("bw.lscvg warns and resets invalid boundary values", {
  x <- seq(0, 2 * pi, length.out = 5)
  expect_cli_warning(
    result <- bw.lscvg(x, lower = -5, upper = 5),
    1,
    "! The boundaries must be positive numbers and 'lower' must be smaller that 'upper'. Default boundaries lower=0, upper=60 were used."
  )
  expect_type(result, "double")

  expect_cli_warning(
    result <- bw.lscvg(x, lower = 10, upper = 5),
    1,
    "! The boundaries must be positive numbers and 'lower' must be smaller that 'upper'. Default boundaries lower=0, upper=60 were used."
  )
  expect_type(result, "double")
})

test_that("bw.lscvg warns when minimum is at edge of the range", {
  x <- rep(0, 10)
  expect_cli_warning(bw.lscvg(x),
                     1,
                     "! Minimum/maximum occurred at one end of the range.")
})

