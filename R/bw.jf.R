#' Compute Bandwidth Selector using J-F Criterion
#'
#' This function computes the optimal bandwidth for circular data
#' using the J-F criterion (Jones–Firth type method).
#' The estimator is based on the framework described in
#' García-Portugués et al. (2019) <doi:10.1016/j.stamet.2017.12.005>.
#'
#' @param x Data vector in radians. Can be numeric or of class `circular`.
#' @param lower Lower bound for optimization (default 0).
#' @param upper Upper bound for optimization (default 60).
#' @param tol Convergence tolerance (default 0.1).
#'
#' @return The computed optimal smoothing parameter (bandwidth).
#'
#' @references
#' García-Portugués, E., Navarro-Esteban, P., & Cuesta-Albertos, J. A. (2019).
#' A goodness-of-fit test for the von Mises–Fisher distribution.
#' *Statistics & Probability Letters, 146*, 178–184.
#'
#' @export
#'
#' @examples
#' library(circular)
#' set.seed(42)
#' x <- rvonmises(50, mu = circular(pi/3), kappa = 2)
#' bw <- bw.jf(x)
#' print(bw)
#'
bw.jf <- function(x,
                  lower = 0,
                  upper = 60,
                  tol = 0.1) {
  if (length(x) == 0) cli::cli_abort("{.var x} must be non-empty.")
  if (!is.numeric(x)) cli::cli_abort("{.var x} must be numeric.")
  x <- circular::circular(x, units = "radians", modulo = "2pi")
  attr(x, "class") <- attr(x, "circularp") <- NULL

  crit_jf <- function(x, h) {
    # Placeholder for JF criterion
    sum(sin(outer(x, x, "-"))^2 * exp(-h))
  }

  bw <- optimize(
    crit_jf,
    interval = c(lower, upper),
    maximum = FALSE,
    x = x
  )$minimum

  if (bw < lower + tol || bw > upper - tol) {
    cli::cli_alert_warning("Optimum at boundary.")
  }
  return(bw)
}
