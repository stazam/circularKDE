#' Adaptive Circular Kernel Density Estimation
#'
#' This function computes an adaptive kernel density estimate for circular data,
#' adjusting the bandwidth locally based on a global bandwidth parameter and data
#' density (see <doi:10.1007/s00180-023-01401-0>). The density is evaluated over
#' a specified range of points between `from` and `to`, using a local adaptation
#' factor to refine the estimate.
#'
#' @param x Data for which the density is to be estimated. The object is coerced to a
#'   numeric vector in radians using `circular::conversion.circular`. Can be a numeric
#'   vector or an object of class `circular`.
#' @param bw0 Global bandwidth parameter, a positive numeric value that sets the
#'   baseline smoothness of the density estimate. Controls the overall scale of the
#'   kernel.
#' @param alpha Numeric scalar between 0 and 1 (default is 0.5) that determines the
#'   sensitivity of the local bandwidth adaptation. Higher values increase the
#'   influence of local data density on the bandwidth.
#' @param type Character string specifying the method for calculating the local
#'   adaptation factor (default is "n"). Options depend on the implementation of the
#'   `local.factor` function (not shown here).
#' @param z Optional numeric vector of points at which to evaluate the density. If
#'   `NULL` (default), a grid of `n` equally spaced points between `from` and `to` is
#'   generated automatically.
#' @param from Starting point of the evaluation range, a numeric value in radians.
#'   Must be finite. Default is `circular(0)` (0 radians).
#' @param to Ending point of the evaluation range, a numeric value in radians. Must
#'   be finite. Default is `circular(2 * pi)` (2π radians).
#' @param n Positive integer specifying the number of evaluation points (default is
#'   500). Ignored if `z` is provided. Determines the resolution of the density
#'   estimate when `z` is `NULL`.
#'
#' @return A numeric vector of length equal to the length of `z` (or `n` if `z` is
#'   `NULL`), containing the estimated density values at each evaluation point.
#'
#' @export
#'
#' @examples
#' # Example with numeric data in radians
#' library(circular)
#' x <- rvonmises(100, mu = circular(0), kappa = 1)
#' bw0 <- bw.lscvg(x = x)
#' dens <- adaptive.density.circular(x, bw0 = bw0)
#' plot(seq(0, 2 * pi, length.out = 500), dens, type = "l",
#'      main = "Adaptive Circular Density")
#'
#' # Example with circular data and custom evaluation points
#' x <- rvonmises(100, mu = circular(0), kappa = 1)
#' bw0 <- bw.lscvg(x = x)
#' z <- seq(0, 2 * pi, length.out = 200)
#' dens <- adaptive.density.circular(x, bw0 = 0.5, z = z)
#' plot(z, dens, type = "l", main = "Density with Custom Points")
#'
#' @import circular
#' @import cli
adaptive.density.circular <- function(x,
                                      bw0,
                                      alpha = 0.5,
                                      type = "n",
                                      z = NULL,
                                      from = 0,
                                      to = 2 * pi,
                                      n = 500) {
  n_x <- length(x)
  if (n_x == 0) {
    cli::cli_abort(
      c("{.var x} must be a non-empty object.", "x" = "You've supplied an object of length {n_x}.")
    )
  }
  if (!is.numeric(x)) {
    if (all(is.na(x))) {
      cli::cli_abort("{.var x} contains all missing values.")
    }
    cli::cli_abort(
      c("{.var x} must be a numeric vector.", "x" = "You've supplied a {.cls {class(x)}} vector.")
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
  attr(x, "class") <- attr(x, "circularp") <- NULL
  if (any(is.na(x))) {
    cli::cli_alert_warning("{.var x} contains missing values, which will be removed.")
    x <- x[!is.na(x)]
  }
  if (!is.numeric(from) | !is.finite(from)) {
    cli_abort("Argument {.var from} must be finite numeric value.")
  }
  if (!is.numeric(to) | !is.finite(to)) {
    cli_abort("Argument {.var to} must be finite numeric value.")
  }
  if (!is.numeric(n) | !is.finite(n)) {
    cli_abort("Argument {.var n} must be finite numeric value.")
  }
  if (round(n) != n | n <= 0) {
    cli_abort("Argument {.var n} must be positive integer.")
  }
  if (is.null(z)) {
    z <- circular(
      seq(
        from = from,
        to = to,
        length.out = n
      ),
      units = "radians",
      zero = 0,
      rotation = "counter",
      modulo = "2pi",
      template = "none"
    )
  }
  kernel.density.adaptive.est <- function(z, x, bw0, alpha, type) {
    n <- length(x)
    factor <- 1 / (2 * n * pi)

    lambda <- local.factor(x, bw0, alpha, type)
    main.part <- sum(1 / besselI(lambda * bw0, 0) * exp (lambda * bw0 * cos(z - x)))
    result <- factor * main.part
    return(result)
  }
  y <- sapply(
    z,
    kernel.density.adaptive.est,
    x = x,
    bw0 = bw0,
    alpha = alpha,
    type = type
  )
  return(y)
}
