test_that("bwCcv returns a numeric value for valid input", {
  set.seed(60)
  x <- rvonmises(50, circular(pi / 2), 1, control.circular = list(units = "radians"))
  result <- bwCcv(x)
  expect_equal(result, 7.98876150)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("bwCcv returns a numeric value for valid input with different seed",
          {
            set.seed(123)
            x <- rvonmises(50, circular(pi / 2), 1)
            result <- bwCcv(x)
            expect_equal(result, 0.808291373)
            expect_type(result, "double")
            expect_length(result, 1)
          })

test_that("bwCcv throws error on empty input", {
  expect_error(bwCcv(numeric(0)), "`x` must be a non-empty object.")
})

test_that("bwCcv throws error if input is not numeric", {
  expect_error(bwCcv(c("a", "b")), "must be a numeric vector")
})

test_that("bwCcv throws error if x contains only NAs", {
  x <- c(NA, NA)
  expect_error(bwCcv(x), "`x` contains all missing values.")
})

test_that("bwCcv removes NA values and returns result", {
  x <- circular(c(0, pi / 2, NA, pi))
  result <- bwCcv(x)
  expect_type(result, "double")
  expect_cli_warning(bwCcv(x),
                     1,
                     "! `x` contains missing values, which will be removed.")
})

test_that("bwCcv handles non-numeric lower", {
  x <- circular(seq(0, 2 * pi, length.out = 5))
  expect_cli_warning(
    result <- bwCcv(x, lower = "zero"),
    1,
    "! Argument `lower` must be numeric. Default value 0 for lower boundary was used."
  )
  expect_type(result, "double")
})

test_that("bwCcv handles non-numeric upper", {
  x <- circular(seq(0, 2 * pi, length.out = 5))
  expect_cli_warning(
    result <- bwCcv(x, upper = "sixty"),
    1,
    "! Argument `upper` must be numeric. Default value 60 for upper boundary was used."
  )
  expect_type(result, "double")
})

test_that("bwCcv warns and resets invalid boundary values", {
  x <- circular(seq(0, 2 * pi, length.out = 5))
  expect_cli_warning(
    result <- bwCcv(x, lower = -5, upper = 5),
    1,
    "! The boundaries must be positive numbers and 'lower' must be smaller than 'upper'. Default boundaries lower=0, upper=60 were used."
  )
  expect_type(result, "double")

  expect_cli_warning(
    result <- bwCcv(x, lower = 10, upper = 5),
    1,
    "! The boundaries must be positive numbers and 'lower' must be smaller than 'upper'. Default boundaries lower=0, upper=60 were used."
  )
  expect_type(result, "double")
})

test_that("bwCcv warns when minimum is at edge of the range", {
  x <- circular(rep(0, 10))
  expect_cli_warning(bwCcv(x, tol = 1),
                     1,
                     "! Minimum/maximum occurred at one end of the range.")
})
