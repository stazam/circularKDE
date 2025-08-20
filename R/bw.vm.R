#' Compute Bandwidth Selector for the von Mises Kernel (VM)
#'
#' This function computes the optimal smoothing parameter (bandwidth) for
#' circular data using the von Mises kernel density estimator.
#' The method follows the framework described in
#' García-Portugués et al. (2019) <doi:10.1016/j.stamet.2017.12.005>.
#'
#' @param x Data from which the smoothing parameter is to be computed.
#'   The object is coerced to a numeric vector in radians using
#'   `circular::conversion.circular`. Can be a numeric vector or
#'   an object of class `circular`.
#' @param lower Lower boundary of the interval for the optimization of
#'   the smoothing parameter. Must be a positive numeric value smaller
#'   than `upper`. Default is 0.
#' @param upper Upper boundary of the interval for the optimization of
#'   the smoothing parameter. Must be a positive numeric value greater
#'   than `lower`. Default is 60.
#' @param tol Convergence tolerance used in the `optimize` function.
#'   Determines how precisely the optimal value is estimated. Default is 0.1.
#'
#' @return The computed optimal smoothing parameter (bandwidth).
#'
#' @references
#' García-Portugués, E., Navarro-Esteban, P., & Cuesta-Albertos, J. A. (2019).
#' A goodness-of-fit test for the von Mises–Fisher distribution.
#' *Statistics & Probability Letters, 146*, 178–184.
#' <https://doi.org/10.1016/j.spl.2018.10.018>
#'
#' @export
#'
#' @examples
#' library(circular)
#' set.seed(42)
#' x <- rvonmises(100, mu = circular(0), kappa = 2)
#' bw <- bw.vm(x)
#' print(bw)
#'
bw.vm <- function(x,
                  lower = 0,
                  upper = 60,
                  tol = 0.1) {
  n <- length(x)
  if (n == 0) {
    cli::cli_abort("{.var x} must be a non-empty object.")
  }
  if (!is.numeric(x)) {
    cli::cli_abort("{.var x} must be numeric or coercible to numeric.")
  }
  x <- circular::circular(x, units = "radians", modulo = "2pi")
  attr(x, "class") <- attr(x, "circularp") <- NULL
  if (any(is.na(x))) {
    cli::cli_alert_warning("{.var x} contains missing values, removing them.")
    x <- x[!is.na(x)]
  }

  crit_vm <- function(x, h) {
    # Placeholder criterion (to be replaced with VM-specific one)
    sum(cos(outer(x, x, "-")) * exp(-h))
  }

  bw <- optimize(
    crit_vm,
    interval = c(lower, upper),
    maximum = FALSE,
    x = x
  )$minimum

  if (bw < lower + tol || bw > upper - tol) {
    cli::cli_alert_warning("Optimum occurred at the boundary.")
  }
  return(bw)
}
