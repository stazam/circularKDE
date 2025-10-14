#' @title Smoothed Cross-Validation for Circular Data 
#'
#' @description This function computes the optimal smoothing parameter (bandwidth) for circular data using a smoothed cross-validation 
#' (SCV) method (see <doi:10.1007/s00180-023-01401-0>). 
#'
#' @param x Data from which the smoothing parameter is to be computed. The object is
#'   coerced to a numeric vector in radians using \code{\link[circular]{conversion.circular}}.
#'   Can be a numeric vector or an object of class \code{circular}.
#' @param np An integer specifying the number of points used in numerical integration
#'   to evaluate the SCV criterion. A higher number increases precision but also
#'   computational cost. Default is 500.
#' @param lower Lower boundary of the interval for the optimization of the smoothing
#'   parameter \code{kappa}. Must be a positive numeric value smaller than \code{upper}.
#'   Default is 0.
#' @param upper Upper boundary of the interval for the optimization of the smoothing
#'   parameter \code{kappa}. Must be a positive numeric value greater than \code{lower}.
#'   Default is 60.
#' @param tol Convergence tolerance for the \code{\link[stats]{optimize}} function, determining the
#'   precision of the optimization process. Also used to detect convergence near boundaries:
#'   if the optimal value is within \code{tol} of \code{lower} or \code{upper}, a warning
#'   is issued suggesting interval adjustment. Default is 0.1.
#' 
#' @details The smoothed cross-validation (SCV) method is an alternative bandwidth 
#' selection approach, originally introduced by Hall & Marron (1992) for linear 
#' densities and adapted for circular data by Zámečník et al. (2023).
#' 
#' The SCV criterion is given by
#' \deqn{\mathrm{SCV}(\kappa) = \frac{R(K)}{nh} 
#'  + \frac{1}{n^{2}} \sum_{i=1}^{n} \sum_{j=1}^{n} 
#'     \big(K_{\kappa} * K_{\kappa} * K_{\kappa} * K_{\kappa} - 2K_{\kappa} * K_{\kappa} *K_{\kappa} + K_{\kappa} * K_{\kappa}\big)(\Theta_i - \Theta_j)}
#' where \eqn{K_\kappa} is the Von Mises kernel with concentration \eqn{\kappa} (for the formula see 3.7, 3.8 in Zámečník et al. (2023)). The optimal bandwidth minimizes the sum 
#' \eqn{ISB(\kappa) + IV(\kappa)} over the interval \code{[lower, upper]}. 
#' 
#' The integral expressions involved in the SCV criterion (see Sections 3.2 in Zámečník et al., 2023) are evaluated numerically using the trapezoidal rule 
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
#' bw <- bwScv(x)
#' print(round(bw$minimum, 2))
#'
#' x <- rvonmises(100, mu = circular(0.5), kappa = 2)
#' bw <- bwScv(x)
#' print(round(bw$minimum, 2))
#'
#' @references
#' Zámečník, S., Horová, I., Katina, S., & Hasilová, K. (2023). An adaptive 
#' method for bandwidth selection in circular kernel density estimation. 
#' \emph{Computational Statistics}.
#' \doi{10.1007/s00180-023-01401-0}
#' 
#' Hall, P., & Marron, J. S. (1992). On the amount of noise inherent in bandwidth 
#' selection for a kernel density estimator. \emph{The Annals of Statistics}, 
#' 20(1), 163-181.
#'
#' @importFrom stats optimize
#' @import circular
#' @import cli
bwScv <- function(x,
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
    modulo = "2pi",
    template = "none"
  )
  x <- as.numeric(conversion.circular(x))
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
    B <- suppressWarnings(besselI(kappa * sqrt(2 * (1 + C)), 0))
    E <- exp(kappa * C)
    b0_kappa <- suppressWarnings(besselI(kappa, 0))

    factor_1 <- 1 / (4 * pi ^ 2 * n * (n - 1) * b0_kappa ^ 4)
    arg_1 <- t(B) %*% B
    part_1 <- factor_1 * h *  (sum(arg_1) - sum(diag(arg_1)))

    factor_2 <- 1 / (2 * n * (n - 1) * pi ^ 2 * b0_kappa ^ 3)
    arg_2 <- t(B) %*% E
    part_2 <- factor_2 * h * (sum(arg_2) - sum(diag(arg_2)))

    factor_3 <- 1 / (2 * n * (n - 1) * pi * b0_kappa ^ 2)
    arg_3 <- suppressWarnings(besselI(kappa * sqrt(2 * (1 + D)), 0))
    part_3 <- factor_3 * (sum(arg_3) - sum(diag(arg_3)))

    part_4 <- suppressWarnings(besselI(2 * kappa, 0)) / (2 * n * pi * b0_kappa ^ 2)

    # ISB + IV
    result <- (part_1 - part_2 + part_3) + part_4
    return(result)
  }
  bw <- optimize(
    scv,
    interval = c(lower, upper),
    maximum = FALSE,
    x = x,
    np = np,
    tol = tol
  )
  if (bw$minimum < lower + tol | bw$minimum > upper - tol) {
    cli::cli_alert_warning("Minimum/maximum occurred at one end of the range.")
  }
  return(bw)
}
