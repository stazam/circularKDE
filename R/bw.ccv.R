#' Compute the Optimal Bandwidth for Circular Data using Complete Cross-Validation
#'
#' This function calculates the optimal smoothing parameter (bandwidth) for circular data
#' using the complete cross-validation (CCV) method (see <doi:10.59170/stattrans-2024-024>). It searches for the value of the
#' smoothing parameter `nu` that minimizes the CCV criterion within a specified interval
#' `[lower, upper]`.
#'
#' @param x Data from which the smoothing parameter is to be computed. The object is
#'   coerced to a numeric vector in radians using `circular::conversion.circular`.
#'   Can be a numeric vector or an object of class `circular`.
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
#'   minimizes the complete cross-validation criterion within the interval
#'   `[lower, upper]`.
#'
#' @export
#'
#' @examples
#' # Example with circular data
#' library(circular)
#' set.seed(123)
#' x <- rwrappednormal(100, mu = circular(2), rho = 0.5)
#' bw <- bw.ccv(x)
#' print(bw)
#'
#' x <- rvonmises(100, mu = circular(0), kappa = 1)
#' bw <- bw.ccv(x)
#' print(bw)
#'
#' @importFrom stats optimize
#' @import circular
#' @import cli
bw.ccv <- function(x,
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
  x <- circular(
    x,
    units = "radians",
    zero = 0,
    rotation = "counter",
    modulo = "2pi"
  )
  attr(x, "class") <- attr(x, "circularp") <- NULL
  if (any(is.na(x))) {
    cli::cli_alert_warning("{.var x} contains missing values, which will be removed")
    x <- x[!is.na(x)]
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

  ccv <- function(x, nu) {
    n <- length(x)
    grid <- outer(x, x, "-")

    b0.nu <- besselI(nu, 0)
    b1.nu <- besselI(nu, 1)
    b2.nu <- besselI(nu, 2)

    cos.grid <- cos(grid)
    sin.grid <- sin(grid)
    exp.nu.cos <- exp(nu * cos.grid)

    cos.0 <- 1
    sin.0 <- 0
    exp.0 <- exp(nu * cos.0)

    factor.1 <- 1 / (2 * pi * n ^ 2 * b0.nu ^ 2)
    part.1 <- factor.1 * sum(besselI(nu * sqrt(2 * (1 + cos.grid)), 0))

    factor.2 <- 1 / (n * (n - 1) * 2 * pi * b0.nu)
    part.2 <- factor.2 * (sum(exp.nu.cos) - n * exp.0)

    factor.3 <- (1 / (2 * nu)) * (b1.nu / b0.nu) * (1 / (2 * pi * b0.nu * n * (n - 1)))
    arg.3 <- exp.nu.cos * (nu ^ 2 * sin.grid ^ 2 - nu * cos.grid)
    arg.3.1 <- exp.0 * (-nu)
    part.3 <- -factor.3 * (sum(arg.3) - n * arg.3.1)

    factor.4 <- (1 / (8 * nu ^ 2)) * (2 * (b1.nu / b0.nu) ^ 2 - b2.nu / b0.nu) * (1 / (2 * pi * b0.nu * n * (n - 1)))
    arg.4 <- exp.nu.cos * (
      nu ^ 4 * sin.grid ^ 4 - 6 * nu ^ 3 * sin.grid ^ 2 * cos.grid +
        3 * nu ^ 2 * (cos.grid ^ 2 - sin.grid ^ 2) - nu ^
        2 * sin.grid ^ 2 + nu * cos.grid
    )
    arg.4.1 <- exp.0 * (3 * nu ^ 2 + nu)
    part.4 <- factor.4 * (sum(arg.4) - n * arg.4.1)

    result <- part.1 - part.2 + part.3 + part.4
    return(result)
  }
  bw <- optimize(
    ccv,
    interval = c(lower, upper),
    maximum = FALSE,
    x = x
  )$minimum
  if (bw < lower + tol | bw > upper - tol) {
    cli::cli_alert_warning("Minimum/maximum occurred at one end of the range.")
  }
  return(bw)
}
