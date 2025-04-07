test_that("bw.scv returns a numeric value for valid input", {
  set.seed(60)
  x <- rvonmises(100, circular(3 * pi / 2), 2, control.circular = list(units = "radians"))
  result <- bw.scv(x)
  expect_equal(result, 6.05465743)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("bw.scv returns a numeric value for valid input with different seed", {
  set.seed(123)
  x <- rvonmises(100, circular(3 * pi / 2), 2, control.circular = list(units = "radians"))
  result <- bw.scv(x)
  expect_equal(result, 4.84740614)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("bw.scv throws error on empty input", {
  expect_error(bw.scv(numeric(0)), "`x` must be a non-empty object.")
})

test_that("bw.scv throws error if input is not numeric", {
  expect_error(bw.scv(c("a", "b")), "must be a numeric vector")
})

test_that("bw.scv throws error if x contains only NAs", {
  x <- c(NA, NA)
  expect_error(bw.scv(x), "`x` contains all missing values.")
})

test_that("bw.scv removes NA values and returns result", {
  x <- c(0, pi / 2, NA, pi)
  result <- bw.scv(x)
  expect_type(result, "double")
  expect_cli_warning(bw.scv(x), 1, "! `x` contains missing values, which will be removed.")
})

test_that("bw.scv handles non-numeric np", {
  x <- seq(0, 2 * pi, length.out = 5)
  expect_cli_warning(
    result <- bw.scv(x, np = "seventy-five"),
    1,
    "! Argument `np` must be numeric. Default value 75 for number of points for evalutaion of numerical integration was used."
  )
  expect_type(result, "double")
})

test_that("bw.scv handles non-numeric lower", {
  x <- seq(0, 2 * pi, length.out = 5)
  expect_cli_warning(
    result <- bw.scv(x, lower = "zero"),
    1,
    "! Argument `lower` must be numeric. Default value 0 for lower boundary was used."
  )
  expect_type(result, "double")
})

test_that("bw.scv handles non-numeric upper", {
  x <- seq(0, 2 * pi, length.out = 5)
  expect_cli_warning(
    result <- bw.scv(x, upper = "sixty"),
    1,
    "! Argument `upper` must be numeric. Default value 60 for upper boundary was used."
  )
  expect_type(result, "double")
})

test_that("bw.scv warns and resets invalid boundary values", {
  x <- seq(0, 2 * pi, length.out = 5)
  expect_cli_warning(
    result <- bw.scv(x, lower = -5, upper = 5),
    1,
    "! The boundaries must be positive numbers and 'lower' must be smaller that 'upper'. Default boundaries lower=0, upper=60 were used."
  )
  expect_type(result, "double")

  expect_cli_warning(
    result <- bw.scv(x, lower = 10, upper = 5),
    1,
    "! The boundaries must be positive numbers and 'lower' must be smaller that 'upper'. Default boundaries lower=0, upper=60 were used."
  )
  expect_type(result, "double")
})

test_that("bw.scv warns when minimum is at edge of the range", {
  x <- rep(0, 10)
  expect_cli_warning(bw.scv(x), 1, "! Minimum/maximum occurred at one end of the range.")
})
