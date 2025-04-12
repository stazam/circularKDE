#' Compute the Optimal Bandwidth for Circular Data using Smoothed Cross-Validation
#'
#' This function computes the optimal smoothing parameter (bandwidth) for circular data
#' using a smoothed cross-validation (SCV) method (see <doi:doi.org/10.1007/s00180-023-01401-0>).
#' It searches for the value of the smoothing parameter `nu` that minimizes the SCV criterion within the
#' specified interval `[lower, upper]`.
#'
#' @param x Data from which the smoothing parameter is to be computed. The object is
#'   coerced to a numeric vector in radians using `circular::conversion.circular`.
#'   Can be a numeric vector or an object of class `circular`.
#' @param np An integer specifying the number of points used in numerical integration
#'   to evaluate the SCV criterion. A higher number increases precision but also
#'   computational cost. Default is 500.
#' @param lower Lower boundary of the interval for the optimization of the smoothing
#'   parameter `mu`. Must be a positive numeric value smaller than `upper`.
#'   Default is 0.
#' @param upper Upper boundary of the interval for the optimization of the smoothing
#'   parameter `mu`. Must be a positive numeric value greater than `lower`.
#'   Default is 60.
#' @param tol Convergence tolerance used in the `optimize` function. Determines how
#'   precisely the optimal value is estimated. Default is 0.1.
#'
#' @return The computed optimal smoothing parameter `mu`, the numeric value
#' that minimizes the smoothed cross-validation criterion.
#'
#' @export
#'
#' @examples
#' # Example with numeric data in radians
#' library(circular)
#  x <- rwrappednormal(100, mu = circular(2), rho = 0.5)
#' bw <- bw.scv(x)
#' print(bw)
#
#' # Example with circular data
#' library(circular)
#' x <- rvonmises(100, mu = circular(0.5), kappa = 2)
#' bw <- bw.scv(x)
#' print(bw)
#'
#' @importFrom circular conversion.circular
#' @importFrom stats optimize
#' @import cli
bw.scv <- function(x,
                   np = 500,
                   lower = 0,
                   upper = 60,
                   tol = 0.1) {
  n <- length(x)
  if (n == 0) {
    cli::cli_abort(
      c("{.var x} must be a non-empty object. ", "x" = "You've supplied an object of length {n}.")
    )
  }
  if (!is.numeric(x)) {
    if (all(is.na(x))) {
      cli::cli_abort("{.var x} contains all missing values.")
    }
    cli::cli_abort(
      c("{.var x} must be a numeric vector. ", "x" = "You've supplied a {.cls {class(x)}} vector.")
    )
  }
  x <- conversion.circular(
    x,
    units = "radians",
    zero = 0,
    rotation = "counter",
    modulo = "2pi"
  )
  attr(x, "class") <- attr(x, "circularp") <- NULL
  if (any(is.na(x))) {
    cli::cli_alert_warning("{.var x} contains missing values, which will be removed.")
    x <- x[!is.na(x)]
  }
  if (!is.numeric(np)) {
    cli::cli_alert_warning(
      c(
        "Argument {.var np} must be numeric. ",
        "Default value 500 for number of points for evalutaion of numerical integration was used."
      )
    )
    np <- 500
  }
  if (!is.numeric(lower)) {
    cli::cli_alert_warning(
      c(
        "Argument {.var lower} must be numeric. ",
        "Default value 0 for lower boundary was used."
      )
    )
    upper <- 0
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
  scv <- function(x, mu, np) {
    # trapezoidal rule for numerical integration in used
    n <- length(x)
    h <- 2 * pi / np
    knots <- seq(0, 2 * pi * (m-1) / m, length = np) # correction for periodic functions

    C <- cos(outer(knots, x, "-"))
    D <- cos(outer(x, x, "-"))
    B <- besselI(mu * sqrt(2 * (1 + C)), 0)
    E <- exp(mu * C)
    bessel_mu <- besselI(mu, 0)

    # first part of ISB integral
    BTB <- t(B) %*% B
    sum_BTB <- sum(BTB) - sum(diag(BTB))
    factor1 <- 1 / (4 * pi ^ 2 * n * (n - 1) * bessel_mu ^ 4)
    first.part <- factor1 * h * sum_BTB

    # second part of ISB integral
    BTE <- t(B) %*% E
    sum_BTE <- sum(BTE) - sum(diag(BTE))
    factor2 <- 1 / (2 * n * (n - 1) * pi ^ 2 * bessel_mu ^ 3)
    second.part <- factor2 * h * sum_BTE

    # third part of ISB integral
    part <- besselI(mu * sqrt(2 * (1 + D)), 0)
    sum_part <- sum(part) - sum(diag(part))
    factor3 <- 1 / (2 * n * (n - 1) * pi * bessel_mu ^ 2)
    third.part <- factor3 * sum_part

    # IV
    last_term <- besselI(2 * mu, 0) / (2 * n * pi * bessel_mu ^ 2)

    # ISB + IV
    scv_value <- first.part - second.part + third.part + last_term
    return(scv_value)
  }
  bw <- optimize(
    scv,
    interval = c(lower, upper),
    maximum = FALSE,
    x = x,
    np = np
  )$minimum
  if (bw < lower + tol | bw > upper - tol) {
    cli::cli_alert_warning("Minimum/maximum occurred at one end of the range.")
  }
  return(bw)
}


