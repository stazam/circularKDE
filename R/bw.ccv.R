#' Compute the Optimal Bandwidth for Circular Data using Complete Cross-Validation
#'
#' This function calculates the optimal smoothing parameter (bandwidth) for circular data
#' using the complete cross-validation (CCV) method (see \doi{10.59170/stattrans-2024-024}). 
#'
#' @param x Data from which the smoothing parameter is to be computed. The object is
#'   coerced to a numeric vector in radians using \code{\link[circular]{conversion.circular}}.
#'   Can be a numeric vector or an object of class \code{circular}.
#' @param lower Lower boundary of the interval to be used in the search for the
#'   smoothing parameter \code{kappa}. Must be a positive numeric value less than \code{upper}.
#'   Default is 0.
#' @param upper Upper boundary of the interval to be used in the search for the
#'   smoothing parameter \code{kappa}. Must be a positive numeric value greater than \code{lower}.
#'   Default is 60.
#' @param tol Convergence tolerance for the \code{\link[stats]{optimize}} function, determining the
#'   precision of the optimization process. Default is 0.1.
#' 
#' @details The complete cross-validation (CCV) method is an alternative for bandwidth 
#' selection, originally proposed by Jones (1991) for linear densities. Its adaptation 
#' to the circular setting was recently studied by Hasilová et al. (2024).
#' 
#' The method uses functionals \eqn{T_m} defined as:
#' \deqn{T_m(\kappa) = \frac{(-1)^m}{n(n-1)}\sum_{i=1}^n\sum_{j \neq i}^n K_{\kappa}^{(2m)}(\theta_{i} - \theta_{j})}
#' where \eqn{K_{\kappa}^{(2m)}} is the (2m)-th derivative of \eqn{K_{\kappa}}.
#' 
#' The CCV criterion can be expressed as:
#' \deqn{CCV(\kappa) = R(f(\kappa)) - T_0(\kappa) + \frac{1}{2}\bar{\sigma}_h^2 T_1(\kappa) + \frac{1}{24}(\eta_{2}^4(K_{\kappa}) - \eta_{4}(K_{\kappa}))T_2(\kappa)}
#' where \eqn{\eta_{j}(K_{\kappa})} denotes the j-th moment of the kernel.
#' 
#' For the von Mises kernel, the explicit CCV criterion becomes:
#' \deqn{CCV(\kappa) = \frac{1}{n^2} \sum_{i=1}^n \sum_{j=1}^n (K_{\kappa} * K_{\kappa})(\theta_i - \theta_j) - T_0(\kappa) + \frac{A_1(\kappa)}{2\kappa}T_1(\kappa) + \frac{2A_1^2(\kappa) - A_2(\kappa)}{8\kappa^2}T_2(\kappa)}
#' where \eqn{A_k(\kappa) = I_k(\kappa)/I_0(\kappa)} is the ratio of modified Bessel functions.
#' 
#' The optimal bandwidth is obtained by minimizing this criterion over the interval 
#' \code{[lower, upper]}.
#'
#' @return The computed optimal smoothing parameter \code{kappa}, a numeric value that
#'   minimizes the complete cross-validation criterion within the interval
#'   \code{[lower, upper]}.
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
#' @references
#' Hasilová, K., Horová, I., Valis, D., & Zámečník, S. (2024).
#' A comprehensive exploration of complete cross-validation for circular data.
#' \emph{Statistics in Transition New Series}, 25(3):1--12. \doi{10.59170/stattrans-2024-024}
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

  ccv <- function(x, kappa) {
    n <- length(x)
    grid <- outer(x, x, "-")

    b0.kappa <- besselI(kappa, 0)
    b1.kappa <- besselI(kappa, 1)
    b2.kappa <- besselI(kappa, 2)

    cos.grid <- cos(grid)
    sin.grid <- sin(grid)
    exp.kappa.cos <- exp(kappa * cos.grid)

    cos.0 <- 1
    sin.0 <- 0
    exp.0 <- exp(kappa * cos.0)

    factor.1 <- 1 / (2 * pi * n ^ 2 * b0.kappa ^ 2)
    part.1 <- factor.1 * sum(besselI(kappa * sqrt(2 * (1 + cos.grid)), 0))

    factor.2 <- 1 / (n * (n - 1) * 2 * pi * b0.kappa)
    part.2 <- factor.2 * (sum(exp.kappa.cos) - n * exp.0)

    factor.3 <- (1 / (2 * kappa)) * (b1.kappa / b0.kappa) * (1 / (2 * pi * b0.kappa * n * (n - 1)))
    arg.3 <- exp.kappa.cos * (kappa ^ 2 * sin.grid ^ 2 - kappa * cos.grid)
    arg.3.1 <- exp.0 * (-kappa)
    part.3 <- -factor.3 * (sum(arg.3) - n * arg.3.1)

    factor.4 <- (1 / (8 * kappa ^ 2)) * (2 * (b1.kappa / b0.kappa) ^ 2 - b2.kappa / b0.kappa) * (1 / (2 * pi * b0.kappa * n * (n - 1)))
    arg.4 <- exp.kappa.cos * (
      kappa ^ 4 * sin.grid ^ 4 - 6 * kappa ^ 3 * sin.grid ^ 2 * cos.grid +
        3 * kappa ^ 2 * (cos.grid ^ 2 - sin.grid ^ 2) - kappa ^
        2 * sin.grid ^ 2 + kappa * cos.grid
    )
    arg.4.1 <- exp.0 * (3 * kappa ^ 2 + kappa)
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
