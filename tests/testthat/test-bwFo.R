test_that("bwFo returns a numeric value for valid input", {
  set.seed(60)
  x <- rvonmises(50, circular(pi / 2), 1, control.circular = list(units = "radians"))
  result <- bwFo(x)
  expect_equal(result, 0.228282328)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("bwFo returns a numeric value for valid input", {
  set.seed(123)
  x <- rvonmises(50, circular(pi / 2), 1)
  result <- bwFo(x)
  expect_equal(result, 0.598002222)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("bwFo throws error on empty input", {
  expect_error(bwFo(numeric(0)), "`x` must be a non-empty object.")
})

test_that("bwFo throws error if input is not numeric", {
  expect_error(bwFo(c("a", "b")), "must be a numeric vector")
})

test_that("bwFo throws error if x contains only NAs", {
  x <- c(NA, NA)
  expect_error(bwFo(x), "`x` contains all missing values.")
})

test_that("bwFo removes NA values and returns error", {
  x <- circular(c(0, pi / 2, NA, pi))
  expect_error(bwFo(x), "`x` must be a numeric vector of length at least 5 after removing missing values.")
})

