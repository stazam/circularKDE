test_that("bwScv returns a numeric value for valid input", {
  set.seed(60)
  x <- rvonmises(100,
                 circular(3 * pi / 2),
                 2,
                 control.circular = list(units = "radians"))
  result <- bwScv(x)
  expect_equal(result, 6.27092315)
  expect_length(result, 1)
  expect_type(result, "double")
})

test_that("bwScv returns a numeric value for valid input with different seed",
          {
            set.seed(123)
            x <- rvonmises(100,
                           circular(3 * pi / 2),
                           2,
                           control.circular = list(units = "radians"))
            result <- bwScv(x)
            expect_equal(result, 5.015259)
            expect_type(result, "double")
            expect_length(result, 1)
          })

test_that("bwScv throws error on empty input", {
  expect_error(bwScv(numeric(0)), "`x` must be a non-empty object.")
})

test_that("bwScv throws error if input is not numeric", {
  expect_error(bwScv(c("a", "b")), "must be a numeric vector")
})

test_that("bwScv throws error if x contains only NAs", {
  x <- c(NA, NA)
  expect_error(bwScv(x), "`x` contains all missing values.")
})

test_that("bwScv removes NA values and returns result", {
  x <- circular(c(0, pi / 2, NA, pi))
  result <- bwScv(x)
  expect_type(result, "double")
  expect_cli_warning(bwScv(x),
                     1,
                     "! `x` contains missing values, which will be removed.")
})

test_that("bwScv handles non-numeric np", {
  x <- circular(seq(0, 2 * pi, length.out = 5))
  expect_cli_warning(
    result <- bwScv(x, np = "seventy-five"),
    1,
    "! Argument `np` must be numeric. Default value 500 for number of points for evaluation of numerical integration was used."
  )
  expect_type(result, "double")
})

test_that("bwScv handles non-numeric lower", {
  x <- circular(seq(0, 2 * pi, length.out = 5))
  expect_cli_warning(
    result <- bwScv(x, lower = "zero"),
    1,
    "! Argument `lower` must be numeric. Default value 0 for lower boundary was used."
  )
  expect_type(result, "double")
})

test_that("bwScv handles non-numeric upper", {
  x <- circular(seq(0, 2 * pi, length.out = 5))
  expect_cli_warning(
    result <- bwScv(x, upper = "sixty"),
    1,
    "! Argument `upper` must be numeric. Default value 60 for upper boundary was used."
  )
  expect_type(result, "double")
})

test_that("bwScv warns and resets invalid boundary values", {
  x <- circular(seq(0, 2 * pi, length.out = 5))
  expect_cli_warning(
    result <- bwScv(x, lower = -5, upper = 5),
    1,
    "! The boundaries must be positive numbers and 'lower' must be smaller than 'upper'. Default boundaries lower=0, upper=60 were used."
  )
  expect_type(result, "double")

  expect_cli_warning(
    result <- bwScv(x, lower = 10, upper = 5),
    1,
    "! The boundaries must be positive numbers and 'lower' must be smaller than 'upper'. Default boundaries lower=0, upper=60 were used."
  )
  expect_type(result, "double")
})

test_that("bwScv warns when minimum is at edge of the range", {
  x <- circular(rep(0, 10))
  expect_cli_warning(bwScv(x),
                     1,
                     "! Minimum/maximum occurred at one end of the range.")
})
