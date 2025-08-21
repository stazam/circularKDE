#' Compute the Optimal Bandwidth for Circular Data using Jones-Faddy Method
#'
#' This function computes the optimal smoothing parameter (bandwidth) for circular data
#' using the Jones-Faddy method. The method is based on higher-order approximations
#' to estimate the optimal bandwidth for kernel density estimation with circular data.
#'
#' @param x Data from which the smoothing parameter is to be computed. The object is
#'   coerced to a numeric vector in radians using `circular::conversion.circular`.
#'   Can be a numeric vector or an object of class `circular`.
#' @param lower Lower boundary of the interval for the optimization of the smoothing
#'   parameter. Must be a positive numeric value smaller than `upper`.
#'   Default is 0.
#' @param upper Upper boundary of the interval for the optimization of the smoothing
#'   parameter. Must be a positive numeric value greater than `lower`.
#'   Default is 60.
#' @param tol Convergence tolerance used in the optimization process. Determines how
#'   precisely the optimal value is estimated. Default is 0.1.
#'
#' @return The computed optimal smoothing parameter, a numeric value derived from
#'   the Jones-Faddy method for circular kernel density estimation.
#'
#' @export
#'
#' @examples
#' # Example with circular data
#' library(circular)
#' set.seed(123)
#' x <- rvonmises(100, mu = circular(0), kappa = 2)
#' bw <- bw.jf(x)
#' print(bw)
#'
#' x <- rwrappednormal(100, mu = circular(1), rho = 0.7)
#' bw <- bw.jf(x)
#' print(bw)
#'
#' @references
#' García-Portugués, E. (2013). Exact risk improvement of bandwidth selectors
#' for kernel density estimation with directional data. \emph{Electronic
#' Journal of Statistics}, 7:1655--1685.
#' \doi{10.1214/13-ejs821}
#'
#' @importFrom stats optimize
#' @import circular
#' @import cli
bw.jf <- function(x,
                  lower = 0,
                  upper = 60,
                  tol = 0.1) {
  if (length(x) == 0) cli::cli_abort("{.var x} must be non-empty.")
  if (!is.numeric(x)) cli::cli_abort("{.var x} must be numeric.")
  x <- circular::circular(x, units = "radians", modulo = "2pi")
  attr(x, "class") <- attr(x, "circularp") <- NULL

  n <- length(x)
  kappa_hat <- circular::mle.vonmises(x)$kappa
  I0 <- besselI(kappa_hat, 0)
  I0k2 <- besselI(2 * kappa_hat, 0)
  I1k2 <- besselI(2 * kappa_hat, 1)
  I2k2 <- besselI(2 * kappa_hat, 2)
  I3k2 <- besselI(2 * kappa_hat, 3)

  R_fVM2 <- (3 * kappa_hat^2 * I2k2 + 2 * kappa_hat * I1k2) / (8 * pi * I0^2)
  R_fVM3 <- (4 * kappa_hat * I1k2 + 30 * kappa_hat^2 * I2k2 + 15 * kappa_hat^3 * I3k2) / (16 * pi * I0^2)
  R_fVM4 <- (8 * kappa_hat^2 * I0k2 + 105 * kappa_hat^4 * I2k2 +
               105 * kappa_hat^3 * I3k2 + 244 * kappa_hat^2 * I2k2) / (32 * pi * I0^2)

  R_hat <- 25 * R_fVM2 / 144 - 5 * R_fVM3 / 36 + R_fVM4 / 36
  kappa_JF <- (16 * sqrt(pi) / 3 * R_hat * n)^(2/9)

  return(kappa_JF)
}
