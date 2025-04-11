test_that("density.adaptive.circular returns a numeric vector for valid input",
          {
            density_vector <- c(
              0.22559, 0.20912, 0.14516, 0.05911, 0.22559
            )
            set.seed(60)
            x <- rvonmises(50, circular(pi / 2), 1, control.circular = list(units = "radians"))
            result <- density.adaptive.circular(x, bw0 = 1, n = 5)
            expect_type(result, "double")
            expect_equal(density_vector, round(result,5))
            expect_length(result, 5)
          })

test_that(
  "density.adaptive.circular returns a numeric vector for valid input with different seed",
  {
    density_vector <- c(
      0.17452, 0.26974, 0.12385, 0.07048, 0.17452
    )
    set.seed(123)
    x <- rvonmises(50, circular(pi / 2), 1)
    result <- density.adaptive.circular(x, bw0 = 1, n = 5)
    expect_equal(density_vector, round(result,5))
    expect_type(result, "double")
    expect_length(result, 5)
  }
)

test_that("density.adaptive.circular throws error on empty input", {
  expect_error(density.adaptive.circular(numeric(0), bw0 = 1),
               "`x` must be a non-empty object.")
})

test_that("density.adaptive.circular throws error if input is not numeric",
          {
            expect_error(density.adaptive.circular(c("a", "b"), bw0 = 1),
                         "must be a numeric vector.")
          })

test_that("density.adaptive.circular throws error if x contains only NAs", {
  x <- c(NA, NA)
  expect_error(density.adaptive.circular(x, bw0 = 1),
               "`x` contains all missing values.")
})

test_that("density.adaptive.circular removes NA values and returns result",
          {
            x <- c(0, pi / 2, NA, pi)
            result <- density.adaptive.circular(x, bw0 = 1)
            expect_type(result, "double")
            expect_length(result, 500)
            expect_cli_warning(
              density.adaptive.circular(x, bw0 = 1),
              1,
              "! `x` contains missing values, which will be removed."
            )
          })

test_that("density.adaptive.circular throws error on non-numeric non-finite from", {
  x <- seq(0, 2 * pi, length.out = 5)
  expect_error(
    density.adaptive.circular(x, bw0 = 1, from = "zero"),
    "Argument `from` must be finite numeric value."
  )
  expect_error(
    density.adaptive.circular(x, bw0 = 1, from = Inf),
    "Argument `from` must be finite numeric value."
  )
})

test_that("density.adaptive.circular throws error on non-numeric to", {
  x <- seq(0, 2 * pi, length.out = 5)
  expect_error(density.adaptive.circular(x, bw0 = 1, to = "two_pi"),
               "Argument `to` must be finite numeric value.")

  expect_error(density.adaptive.circular(x, bw0 = 1, to = NA),
               "Argument `to` must be finite numeric value.")
  }
)


# This one
test_that("density.adaptive.circular throws error on invalid n", {
  x <- seq(0, 2 * pi, length.out = 5)
  expect_error(
    density.adaptive.circular(x, bw0 = 1, n = "five_hundred"),
    "Argument `n` must be finite numeric value."
  )
  expect_error(
    density.adaptive.circular(x, bw0 = 1, n = Inf),
    "Argument `n` must be finite numeric value."
  )
  expect_error(
    density.adaptive.circular(x, bw0 = 1, n = -1),
    "Argument `n` must be positive integer."
  )
  expect_error(
    density.adaptive.circular(x, bw0 = 1, n = 3.12),
    "Argument `n` must be positive integer."
  )
})

test_that("density.adaptive.circular uses custom z correctly", {
  x <- seq(0, 2 * pi, length.out = 5)
  z <- seq(0, pi, length.out = 10)
  result <- density.adaptive.circular(x, bw0 = 1, z = z)
  expect_type(result, "double")
  expect_length(result, 10)
  expect_true(all(is.finite(result)))
})
