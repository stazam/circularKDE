test_that("density.adaptive.circular returns a numeric vector for valid input",
          {
            density_vector <- c(
              0.22558994, 0.24065143, 0.25126071, 0.25704007, 0.25814572,
              0.25521965, 0.24926659, 0.24148486, 0.23308401, 0.22511913,
              0.21836329, 0.21323103, 0.20975753, 0.20763272, 0.20628480,
              0.20500302, 0.20308266, 0.19996752, 0.19535936, 0.18926535,
              0.18196666, 0.17391298, 0.16557292, 0.15728759, 0.14917723,
              0.14113349, 0.13290015, 0.12421165, 0.11493727, 0.10517665,
              0.09527289, 0.08574149, 0.07714614, 0.06996959, 0.06452661,
              0.06094593, 0.05922033, 0.05929951, 0.06118774, 0.06501009,
              0.07102363, 0.07956768, 0.09096334, 0.10538569, 0.12273937,
              0.14256945, 0.16403438, 0.18595609, 0.20694629, 0.22558994
            )
            set.seed(60)
            x <- rvonmises(50, circular(pi / 2), 1, control.circular = list(units = "radians"))
            result <- density.adaptive.circular(x, bw0 = 1, n = 50)
            expect_type(result, "double")
            expect_equal(density_vector, round(result,8))
            expect_length(result, 50)
          })

test_that(
  "density.adaptive.circular returns a numeric vector for valid input with different seed",
  {
    density_vector <- c(
      0.17452359, 0.19060789, 0.20696619, 0.22293133, 0.23781630,
      0.25099414, 0.26195674, 0.27034603, 0.27596072, 0.27874575,
      0.27877165, 0.27620863, 0.27129844, 0.26432665, 0.25559827,
      0.24541971, 0.23408843, 0.22188857, 0.20908794, 0.19593080,
      0.18262434, 0.16932310, 0.15612175, 0.14306726, 0.13019514,
      0.11758224, 0.10539787, 0.09393035, 0.08357219, 0.07476087,
      0.06788890, 0.06320999, 0.06077140, 0.06039595, 0.06172159,
      0.06428842, 0.06764765, 0.07146006, 0.07555502, 0.07993524,
      0.08473342, 0.09014560, 0.09637290, 0.10359311, 0.11196188,
      0.12162288, 0.13270117, 0.14526684, 0.15927729, 0.17452359
    )
    set.seed(123)
    x <- rvonmises(50, circular(pi / 2), 1)
    result <- density.adaptive.circular(x, bw0 = 1, n = 50)
    expect_equal(density_vector, round(result,8))
    expect_type(result, "double")
    expect_length(result, 50)
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
