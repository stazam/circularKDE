test_that("adaptive.density.circular returns a numeric vector for valid input",
          {
            density_vector <- c(0.212165846,
                                0.211115874,
                                0.140838460,
                                0.072948386,
                                0.212165846)
            set.seed(60)
            x <- rvonmises(50, circular(pi / 2), 1)
            result <- adaptive.density.circular(x, bw0 = 1, n = 5)
            expect_type(result, "double")
            expect_equal(density_vector, result)
            expect_length(result, 5)
          })

test_that(
  "adaptive.density.circular returns a numeric vector for valid input with different seed",
  {
    density_vector <- c(0.17507659, 0.26227211, 0.12427576, 0.07484904, 0.17507659)
    set.seed(123)
    x <- rvonmises(50, circular(pi / 2), 1)
    result <- adaptive.density.circular(x, bw0 = 1, n = 5)
    expect_equal(density_vector, result)
    expect_type(result, "double")
    expect_length(result, 5)
  }
)

test_that("adaptive.density.circular throws error on empty input", {
  expect_error(adaptive.density.circular(numeric(0), bw0 = 1),
               "`x` must be a non-empty object.")
})

test_that("adaptive.density.circular throws error if input is not numeric",
          {
            expect_error(adaptive.density.circular(c("a", "b"), bw0 = 1),
                         "must be a numeric vector.")
          })

test_that("adaptive.density.circular throws error if x contains only NAs", {
  x <- circular(c(NA, NA))
  expect_error(adaptive.density.circular(x, bw0 = 1),
               "`x` contains all missing values.")
})

test_that("adaptive.density.circular removes NA values and returns result",
          {
            x <- circular(c(0, pi / 2, NA, pi))
            result <- adaptive.density.circular(x, bw0 = 1)
            expect_type(result, "double")
            expect_length(result, 500)
            expect_cli_warning(
              adaptive.density.circular(x, bw0 = 1),
              1,
              "! `x` contains missing values, which will be removed."
            )
          })

test_that("adaptive.density.circular throws error on non-numeric non-finite from",
          {
            x <- circular(seq(0, 2 * pi, length.out = 5))
            expect_error(
              adaptive.density.circular(x, bw0 = 1, from = "zero"),
              "Argument `from` must be finite numeric value."
            )
            expect_error(
              adaptive.density.circular(x, bw0 = 1, from = Inf),
              "Argument `from` must be finite numeric value."
            )
          })

test_that("adaptive.density.circular throws error on non-numeric to", {
  x <- circular(seq(0, 2 * pi, length.out = 5))
  expect_error(
    adaptive.density.circular(x, bw0 = 1, to = "two_pi"),
    "Argument `to` must be finite numeric value."
  )

  expect_error(
    adaptive.density.circular(x, bw0 = 1, to = NA),
    "Argument `to` must be finite numeric value."
  )
})

test_that("adaptive.density.circular throws error on invalid n", {
  x <- circular(seq(0, 2 * pi, length.out = 5))
  expect_error(
    adaptive.density.circular(x, bw0 = 1, n = "five_hundred"),
    "Argument `n` must be finite numeric value."
  )
  expect_error(
    adaptive.density.circular(x, bw0 = 1, n = Inf),
    "Argument `n` must be finite numeric value."
  )
  expect_error(adaptive.density.circular(x, bw0 = 1, n = -1),
               "Argument `n` must be positive integer.")
  expect_error(
    adaptive.density.circular(x, bw0 = 1, n = 3.12),
    "Argument `n` must be positive integer."
  )
})

test_that("adaptive.density.circular uses custom z correctly", {
  x <- circular(seq(0, 2 * pi, length.out = 5))
  z <- circular(seq(0, pi, length.out = 10))
  result <- adaptive.density.circular(x, bw0 = 1, z = z)
  expect_type(result, "double")
  expect_length(result, 10)
  expect_true(all(is.finite(result)))
})

