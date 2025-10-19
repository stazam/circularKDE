#' @title Adaptive Kernel Density Estimation for circular data
#'
#' @description This function computes an adaptive kernel density estimate for circular data,
#' adjusting the bandwidth locally based on a global bandwidth parameter and data
#' density (see \doi{10.1007/s00180-023-01401-0}).
#'
#' @param x Data for which the density is to be estimated. The object is coerced to a
#'   numeric vector in radians using \code{\link[circular]{circular}}. Can be a numeric
#'   vector or an object of class \code{circular}.
#' @param bw0 Global bandwidth parameter, a positive numeric value that sets the
#'   baseline smoothness of the density estimate. Controls the overall scale of the
#'   kernel.
#' @param alpha Numeric scalar between 0 and 1 (default is 0.5) that determines the
#'   sensitivity of the local bandwidth adaptation. Higher values increase the
#'   influence of local data density on the bandwidth.
#' @param type Character string specifying the method for calculating the local
#'   adaptation factor (default is "n"). Available options are:
#'   \describe{
#'     \item{"am"}{Arithmetic mean of pilot density estimates}
#'     \item{"gm"}{Geometric mean of pilot density estimates (requires all positive values)}
#'     \item{"rv"}{Range (max - min) of pilot density estimates}
#'     \item{"n"}{No scaling factor (equivalent to fixed bandwidth)}
#'   }
#' @param z Optional numeric vector of points at which to evaluate the density. If
#'   \code{NULL} (default), a grid of \code{n} equally spaced points between \code{from} and \code{to} is
#'   generated automatically.
#' @param from Starting point of the evaluation range, a numeric value in radians.
#'   Must be finite. Default is \code{circular(0)} (0 radians).
#' @param to Ending point of the evaluation range, a numeric value in radians. Must
#'   be finite. Default is \code{circular(2 * pi)} (2*pi radians).
#' @param n Positive integer specifying the number of evaluation points (default is
#'   500). Ignored if \code{z} is provided. Determines the resolution of the density
#'   estimate when \code{z} is \code{NULL}.
#'
#' @details The method extends the classical idea of variable bandwidth estimators from the linear case (Breiman et al., 1977) to the circular setting by employing the von Mises kernel. Specifically, the procedure follows three steps:
#' \enumerate{
#' \item A pilot density estimate is obtained using a fixed global bandwidth \eqn{\kappa} (\code{bw0}).
#' \item Local adaptation factors \eqn{\lambda_i} are computed at each data point \eqn{\Theta_i} according to
#' \deqn{\lambda_i = \left\{ g / \hat{f}(\Theta_i) \right\}^{\alpha},}
#' where \eqn{g} is a global measure of central tendency (typically the geometric mean) and \eqn{\alpha \in [0,1]} is a sensitivity parameter.
#' \item The adaptive density is evaluated using local bandwidths \eqn{\lambda_i \cdot \kappa}, resulting in the estimator
#' \deqn{\hat{f}(\theta) = \frac{1}{2n\pi} \sum_{i=1}^{n} \frac{1}{I_0(\lambda_i \kappa)}
#' \exp\!\big(\lambda_i \kappa \cos(\theta - \Theta_i)\big).}
#' }
#'
#' @return A numeric vector of length equal to the length of \code{z} (or \code{n} if \code{z} is
#'   \code{NULL}), containing the estimated density values at each evaluation point.
#'
#' @export
#'
#' @examples
#' # Example with numeric data in radians
#' library(circular)
#' x <- rvonmises(100, mu = circular(0), kappa = 1)
#' bw0 <- bwLscvg(x = x)
#' dens <- adaptiveDensityCircular(x, bw0 = bw0)
#' plot(seq(0, 2 * pi, length.out = 500), dens, type = "l",
#'      main = "Adaptive Circular Density")
#'
#' # Example with numerical integration over interval [0,2\pi] to verify normalization of computed density
#' library(circular)
#' x <- rvonmises(100, mu = circular(0), kappa = 1)
#' bw0 <- bwLscvg(x = x)
#' dens <- function(z)adaptiveDensityCircular(z = z, bw0 = bw0, x=x)
#' integrate(dens, lower = 0, upper = 2*pi) # 1 with absolute error < 4e-07
#'
#' @references
#' Zámečník, S., Horová, I., Katina, S., & Hasilová, K. (2023). An adaptive
#' method for bandwidth selection in circular kernel density estimation.
#' \emph{Computational Statistics}.
#' \doi{10.1007/s00180-023-01401-0}
#'
#' Breiman, L., Meisel, W., & Purcell, E. (1977). Variable kernel estimates of
#' multivariate densities. \emph{Technometrics}, 19(2), 135-144.
#' \doi{10.2307/1268623}
#'
#' @importFrom circular circular
#' @import cli
adaptiveDensityCircular <- function(x,
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
  x <- as.numeric(x)
  if (any(is.na(x))) {
    cli::cli_alert_warning("{.var x} contains missing values, which will be removed.")
    x <- x[!is.na(x)]
  }
  if (bw0 <= 0 | !is.numeric(bw0)) {
    cli::cli_abort("Argument {.var bw0} must be a positive numeric value.")
  }
  if (!is.numeric(alpha) | alpha < 0 | alpha > 1) {
    cli::cli_abort("Argument {.var alpha} must be a numeric value between 0 and 1.")
  }
  if (!is.numeric(from) | !is.finite(from)) {
    cli::cli_abort("Argument {.var from} must be finite numeric value.")
  }
  if (!is.numeric(to) | !is.finite(to)) {
    cli::cli_abort("Argument {.var to} must be finite numeric value.")
  }
  if (!is.numeric(n) | !is.finite(n)) {
    cli::cli_abort("Argument {.var n} must be finite numeric value.")
  }
  if (round(n) != n | n <= 0) {
    cli::cli_abort("Argument {.var n} must be positive integer.")
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
    z <- as.numeric(z)
  }
  # compute the local adaptation factors at every point in x
  lambda <- localFactor(x, bw0, alpha, type)
  kernelDensityAdaptiveEst <- function(x, z, bw0, alpha, type, lambda) {
    n <- length(x)
    factor <- 1 / (2 * n * pi)
    main_part <- sum(1 / besselI(lambda * bw0, 0) * exp (lambda * bw0 * cos(z - x)), na.rm = TRUE)
    result <- factor * main_part
    return(result)
  }
  y <- sapply(
    X = z,
    FUN = kernelDensityAdaptiveEst,
    x = x,
    bw0 = bw0,
    alpha = alpha,
    type = type,
    lambda = lambda
  )
  return(y)
}


