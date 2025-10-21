test_that("adaptiveDensityCircular returns a numeric vector for valid input",
          {
            density_vector <- c(0.225589943,
                                0.209118581,
                                0.145161036,
                                0.059112002,
                                0.225589943)
            set.seed(60)
            x <- rvonmises(50, circular(pi / 2), 1)
            result <- adaptiveDensityCircular(x, bw0 = 1, n = 5)
            expect_type(result, "double")
            expect_equal(density_vector, result)
            expect_length(result, 5)
          })

test_that(
  "adaptiveDensityCircular returns a numeric vector for valid input with different seed",
  {
    density_vector <- c(0.17452359, 0.26973702, 0.12384836, 0.07047693, 0.17452359)
    set.seed(123)
    x <- rvonmises(50, circular(pi / 2), 1)
    result <- adaptiveDensityCircular(x, bw0 = 1, n = 5)
    expect_equal(density_vector, result)
    expect_type(result, "double")
    expect_length(result, 5)
  }
)

test_that("adaptiveDensityCircular throws error on empty input", {
  expect_error(adaptiveDensityCircular(numeric(0), bw0 = 1),
               "`x` must be a non-empty object.")
})

test_that("adaptiveDensityCircular throws error if input is not numeric",
          {
            expect_error(adaptiveDensityCircular(c("a", "b"), bw0 = 1),
                         "must be a numeric vector.")
          })

test_that("adaptiveDensityCircular throws error if x contains only NAs", {
  x <- circular(c(NA, NA))
  expect_error(adaptiveDensityCircular(x, bw0 = 1),
               "`x` contains all missing values.")
})

test_that("adaptiveDensityCircular removes NA values and returns result",
          {
            x <- circular(c(0, pi / 2, NA, pi))
            result <- adaptiveDensityCircular(x, bw0 = 1)
            expect_type(result, "double")
            expect_length(result, 500)
            expect_cli_warning(
              adaptiveDensityCircular(x, bw0 = 1),
              1,
              "! `x` contains missing values, which will be removed."
            )
          })

test_that("adaptiveDensityCircular throws error on non-numeric non-finite from",
          {
            x <- circular(seq(0, 2 * pi, length.out = 5))
            expect_error(
              adaptiveDensityCircular(x, bw0 = 1, from = "zero"),
              "Argument `from` must be finite numeric value."
            )
            expect_error(
              adaptiveDensityCircular(x, bw0 = 1, from = Inf),
              "Argument `from` must be finite numeric value."
            )
          })

test_that("adaptiveDensityCircular throws error on non-numeric to", {
  x <- circular(seq(0, 2 * pi, length.out = 5))
  expect_error(
    adaptiveDensityCircular(x, bw0 = 1, to = "two_pi"),
    "Argument `to` must be finite numeric value."
  )

  expect_error(
    adaptiveDensityCircular(x, bw0 = 1, to = NA),
    "Argument `to` must be finite numeric value."
  )
})

test_that("adaptiveDensityCircular throws error on invalid n", {
  x <- circular(seq(0, 2 * pi, length.out = 5))
  expect_error(
    adaptiveDensityCircular(x, bw0 = 1, n = "five_hundred"),
    "Argument `n` must be finite numeric value."
  )
  expect_error(
    adaptiveDensityCircular(x, bw0 = 1, n = Inf),
    "Argument `n` must be finite numeric value."
  )
  expect_error(adaptiveDensityCircular(x, bw0 = 1, n = -1),
               "Argument `n` must be positive integer.")
  expect_error(
    adaptiveDensityCircular(x, bw0 = 1, n = 3.12),
    "Argument `n` must be positive integer."
  )
})

test_that("adaptiveDensityCircular uses custom z correctly", {
  x <- circular(seq(0, 2 * pi, length.out = 5))
  z <- circular(seq(0, pi, length.out = 10))
  result <- adaptiveDensityCircular(x, bw0 = 1, z = z)
  expect_type(result, "double")
  expect_length(result, 10)
  expect_true(all(is.finite(result)))
})

