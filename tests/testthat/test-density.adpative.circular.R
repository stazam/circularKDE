test_that("density.adaptive.circular returns a numeric vector for valid input",
          {
            set.seed(60)
            x <- rvonmises(50, circular(pi / 2), 1, control.circular = list(units = "radians"))
            result <- density.adaptive.circular(x, bw0 = 1)
            expect_type(result, "double")
            expect_length(result, 500)
            expect_true(all(is.finite(result)))
          })

test_that(
  "density.adaptive.circular returns a numeric vector for valid input with different seed",
  {
    set.seed(123)
    x <- rvonmises(50, circular(pi / 2), 1)
    result <- density.adaptive.circular(x, bw0 = 1)
    expect_type(result, "double")
    expect_length(result, 500)
    expect_true(all(is.finite(result)))
  }
)

test_that("density.adaptive.circular throws error on empty input", {
  expect_error(density.adaptive.circular(numeric(0), bw0 = 1),
               "`x` must be a non-empty object.")
})

test_that("density.adaptive.circular throws error if input is not numeric",
          {
            expect_error(density.adaptive.circular(c("a", "b"), bw0 = 1),
                         "must be a numeric vector")
          })

test_that("density.adaptive.circular throws error if x contains only NAs", {
  x <- c(NA, NA)
  expect_error(density.adaptive.circular(x, bw0 = 1),
               "`x` is after removal of length")
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
              "! `x` contains missing values, which will be removed"
            )
          })

test_that("density.adaptive.circular throws error on non-numeric from", {
  x <- seq(0, 2 * pi, length.out = 5)
  expect_error(
    density.adaptive.circular(x, bw0 = 1, from = "zero"),
    "`from` must be a numeric argument"
  )
})

test_that("density.adaptive.circular throws error on non-finite from", {
  x <- seq(0, 2 * pi, length.out = 5)
  expect_error(
    density.adaptive.circular(x, bw0 = 1, from = Inf),
    "`from` must be a non-finite argument"
  )
})

test_that("density.adaptive.circular throws error on non-numeric to", {
  x <- seq(0, 2 * pi, length.out = 5)
  expect_error(density.adaptive.circular(x, bw0 = 1, to = "two_pi"),
               "`to` must be a numeric argument")
})

test_that("density.adaptive.circular throws error on non-finite to", {
  x <- seq(0, 2 * pi, length.out = 5)
  expect_error(density.adaptive.circular(x, bw0 = 1, to = NA),
               "`to` must be a non-finite argument")
})

test_that("density.adaptive.circular throws error on invalid n", {
  x <- seq(0, 2 * pi, length.out = 5)
  expect_error(
    density.adaptive.circular(x, bw0 = 1, n = "five_hundred"),
    "`n` must be a numeric argument"
  )
  expect_error(
    density.adaptive.circular(x, bw0 = 1, n = -1),
    "argument 'n' must be integer and positive"
  )
})

test_that("density.adaptive.circular uses custom z correctly", {
  x <- seq(0, 2 * pi, length.out = 5)
  z <- seq(0, pi, length.out = 10)
  result <- density.adaptive.circular(x, bw0 = 1, z = z)
  expect_type(result, "double")
  expect_length(result, 10)  # Matches length of z
  expect_true(all(is.finite(result)))
})
