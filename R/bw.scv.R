#' @title Compute the Optimal Bandwidth for Circular Data using Smoothed Cross-Validation
#'
#' @description This function computes the optimal smoothing parameter (bandwidth) for circular data using a smoothed cross-validation 
#' (SCV) method (see \doi{10.1007/s00180-023-01401-0}). 
#'
#' @param x Data from which the smoothing parameter is to be computed. The object is
#'   coerced to a numeric vector in radians using `circular::conversion.circular`.
#'   Can be a numeric vector or an object of class `circular`.
#' @param np An integer specifying the number of points used in numerical integration
#'   to evaluate the SCV criterion. A higher number increases precision but also
#'   computational cost. Default is 500.
#' @param lower Lower boundary of the interval for the optimization of the smoothing
#'   parameter `kappa`. Must be a positive numeric value smaller than `upper`.
#'   Default is 0.
#' @param upper Upper boundary of the interval for the optimization of the smoothing
#'   parameter `kappa`. Must be a positive numeric value greater than `lower`.
#'   Default is 60.
#' @param tol Convergence tolerance used in the `optimize` function. Determines how
#'   precisely the optimal value is estimated. Default is 0.1.
#' 
#' @details The SCV criterion is given by
#' \deqn{\mathrm{SCV}(\kappa) = \frac{R(K)}{nh} 
#'  + \frac{1}{n^{2}} \sum_{i=1}^{n} \sum_{j=1}^{n} 
#'     \big(K_{\kappa} * K_{\kappa} * K_{\kappa} * K_{\kappa} - 2K_{\kappa} * K_{\kappa} *K_{\kappa} + K_{\kappa} * K_{\kappa}\big)(\Theta_i - \Theta_j)}
#' where \eqn{K_\kappa} is the Von Mises kernel with concentration \eqn{\kappa}. The function searches for the value of \eqn{\kappa} 
#' that minimizes this criterion over the interval \code{[lower, upper]}. 
#' 
#' The integral expressions involved in the SCV criterion (see 3.8 in <doi:10.1007/s00180-023-01401-0>) are evaluated numerically using the trapezoidal rule 
#' on a uniform grid of length \code{np}. 
#'
#' @return The computed optimal smoothing parameter \code{kappa}, the numeric value
#' that minimizes the smoothed cross-validation criterion.
#'
#' @export
#'
#' @examples
#' # Example with circular data
#' library(circular)
#' x <- rwrappednormal(100, mu = circular(2), rho = 0.5)
#' bw <- bw.scv(x)
#' print(bw)
#
#' x <- rvonmises(100, mu = circular(0.5), kappa = 2)
#' bw <- bw.scv(x)
#' print(bw)
#'
#' @references
#' Zámečník, S., Horová, I., Katina, S., & Hasilová, K. (2023). An adaptive 
#' method for bandwidth selection in circular kernel density estimation. 
#' \emph{Computational Statistics}.
#' \doi{10.1007/s00180-023-01401-0}
#'
#' @importFrom stats optimize
#' @import circular
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
  x <- circular(
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
        "Default value 500 for number of points for evaluation of numerical integration was used."
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
        "The boundaries must be positive numbers and 'lower' must be smaller than 'upper'. ",
        "Default boundaries lower=0, upper=60 were used."
      )
    )
    lower <- 0
    upper <- 60
  }
  scv <- function(x, kappa, np) {
    # trapezoidal rule for numerical integration is used
    n <- length(x)
    h <- 2 * pi / np
    knots <- seq(0, 2 * pi * (np - 1) / np, length = np) # correction for periodic functions

    C <- cos(outer(knots, x, "-"))
    D <- cos(outer(x, x, "-"))
    B <- besselI(kappa * sqrt(2 * (1 + C)), 0)
    E <- exp(kappa * C)
    b0.kappa <- besselI(kappa, 0)

    factor.1 <- 1 / (4 * pi ^ 2 * n * (n - 1) * b0.kappa ^ 4)
    arg.1 <- t(B) %*% B
    part.1 <- factor.1 * h *  (sum(arg.1) - sum(diag(arg.1)))

    factor.2 <- 1 / (2 * n * (n - 1) * pi ^ 2 * b0.kappa ^ 3)
    arg.2 <- t(B) %*% E
    part.2 <- factor.2 * h * (sum(arg.2) - sum(diag(arg.2)))

    factor.3 <- 1 / (2 * n * (n - 1) * pi * b0.kappa ^ 2)
    arg.3 <- besselI(kappa * sqrt(2 * (1 + D)), 0)
    part.3 <- factor.3 * (sum(arg.3) - sum(diag(arg.3)))

    part.4 <- besselI(2 * kappa, 0) / (2 * n * pi * b0.kappa ^ 2)

    # ISB + IV
    result <- (part.1 - part.2 + part.3) + part.4
    return(result)
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
