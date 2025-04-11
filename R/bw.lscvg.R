#' Compute the Optimal Bandwidth for Circular Data using Generalized Least Squares Cross-Validation
#'
#' This function computes the optimal smoothing parameter (bandwidth) for circular data
#' using a generalized least squares cross-validation method (see <doi:10.1007/s00180-023-01401-0>).
#' It searches for the value of the smoothing parameter `nu` that minimizes the cross-validation criterion within the
#' specified interval `[lower, upper]`.
#'
#' @param x Data from which the smoothing parameter is to be computed. The object is
#'   coerced to a numeric vector in radians using `circular::conversion.circular`.
#'   Can be a numeric vector or an object of class `circular`.
#' @param g A numeric scalar that controls the variability in the cross-validation
#'   procedure. It influences the scaling in the internal calculations, affecting the
#'   bandwidth estimation. Default is 4.
#' @param lower Lower boundary of the interval to be used in the search for the
#'   smoothing parameter `nu`. Must be a positive numeric value less than `upper`.
#'   Default is 0.
#' @param upper Upper boundary of the interval to be used in the search for the
#'   smoothing parameter `nu`. Must be a positive numeric value greater than `lower`.
#'   Default is 60.
#' @param tol Convergence tolerance for the `optimize` function, determining the
#'   precision of the optimization process. Default is 0.1.
#'
#' @return The computed optimal smoothing parameter `nu`, a numeric value that
#'   minimizes the least squares cross-validation criterion within the interval
#'   `[lower, upper]`.
#'
#' @export
#'
#' @examples
#' # Example with numeric data in radians
#' set.seed(123)
#' x <- runif(100, 0, 2 * pi)
#' bw <- bw.lscvg(x)
#' print(bw)
#'
#' # Example with circular data
#' library(circular)
#' x_circ <- rvonmises(100, mu = circular(0), kappa = 1)
#' bw <- bw.lscvg(x_circ)
#' print(bw)
#'
#' @importFrom circular conversion.circular
#' @importFrom stats optimize
#' @import cli
bw.lscvg <- function(x,
                     g = 4,
                     lower = 0,
                     upper = 60,
                     tol = 0.1) {
  n <- length(x)
  if (n == 0) {
    cli::cli_abort(
      c("{.var x} must be a non-empty object.", "x" = "You've supplied an object of length {n}.")
    )
  }
  if (!is.numeric(x)) {
    if (all(is.na(x))) {
      cli::cli_abort("{.var x} contains all missing values.")
    }
    cli::cli_abort(
      c("{.var x} must be a numeric vector", "x" = "You've supplied a {.cls {class(x)}} vector.")
    )
  }
  x <- conversion.circular(
    x,
    units = "radians",
    zero = 0,
    rotation = "counter",
    modulo = "2pi",
    template = "none"
  )
  attr(x, "class") <- attr(x, "circularp") <- NULL
  if (any(is.na(x))) {
    cli::cli_alert_warning("{.var x} contains missing values, which will be removed")
    x <- x[!is.na(x)]
  }
  if (!is.numeric(g)) {
    cli::cli_alert_warning(c(
      "Argument {.var g} must be numeric. ",
      "Default value 4 for coefficient was used."
    ))
    g <- 4
  }
  if (!is.numeric(lower)) {
    cli::cli_alert_warning(
      c(
        "Argument {.var lower} must be numeric. ",
        "Default value 0 for lower boundary was used."
      )
    )
    lower <- 0
  }
  if (!is.numeric(upper)) {
    cli::cli_alert_warning(
      c(
        "Argument {.var upper} must be numeric. ",
        "Default value 60 for upper boundary was used."
      )
    )
    upper <- 60
  }
  if (lower < 0 | lower >= upper) {
    cli::cli_alert_warning(
      c(
        "The boundaries must be positive numbers and 'lower' must be smaller that 'upper'. ",
        "Default boundaries lower=0, upper=60 were used."
      )
    )
    lower <- 0
    upper <- 60
  }

  lscvg <- function(x, nu, g) {
    n <- length(x)
    grid <- outer(x, x, '-')
    cos_grid <- cos(grid)

    b0_nu      <- besselI(nu, 0)
    b0_2nu     <- besselI(2 * nu, 0)
    b0_nu_g    <- besselI(nu / g, 0)
    b0_nu_half <- besselI(nu / 2, 0)

    part1 <- b0_2nu / (2 * pi * n * (b0_nu^2))

    s1 <- sum(exp((nu/g) * cos_grid)) - n * exp(nu / g)
    part2 <- (2 / (g * (g - 2))) * (1 / (2 * pi * b0_nu_g)) * s1

    sqrt_term <- nu * sqrt(2 * (1 + cos_grid))
    s2 <- sum(besselI(sqrt_term, 0)) - n * b0_2nu
    part3 <- ((n - 1) / n) * (1 / (2 * pi * (b0_nu^2))) * s2

    s3 <- sum(exp((nu/2) * cos_grid)) - n * exp(nu / 2)
    part4 <- ((g - 1) / (g - 2)) * (1 / (2 * pi * b0_nu_half)) * s3

    output <- part1 + (part2 + part3 - part4) / (n * (n - 1))
    return(output)
  }
  bw <- optimize(
    lscvg,
    interval = c(lower, upper),
    maximum = FALSE,
    x = x ,
    g = g
  )$minimum
  if (bw < lower + tol | bw > upper - tol) {
    cli::cli_alert_warning("Minimum/maximum occurred at one end of the range.")
  }
  return(bw)
}
