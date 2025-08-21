#' Compute the Optimal Bandwidth for Circular Data using Jones-Faddy Method
#'
#' This function computes the optimal smoothing parameter (bandwidth) for circular data
#' using the Jones-Faddy method. The method is based on higher-order approximations
#' to estimate the optimal bandwidth for kernel density estimation with circular data.
#'
#' @param x Data from which the smoothing parameter is to be computed. The object is
#'   coerced to a numeric vector in radians using `circular::conversion.circular`.
#'   Can be a numeric vector or an object of class `circular`.
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
      c("{.var x} must be a numeric vector.", "x" = "You've supplied a {.cls {class(x)}} vector.")
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

  n <- length(x)
  kappa.hat <- circular::mle.vonmises(x)$kappa
  I0 <- besselI(kappa.hat, 0)
  I0k2 <- besselI(2 * kappa.hat, 0)
  I1k2 <- besselI(2 * kappa.hat, 1)
  I2k2 <- besselI(2 * kappa.hat, 2)
  I3k2 <- besselI(2 * kappa.hat, 3)

  R.fVM2 <- (3 * kappa.hat^2 * I2k2 + 2 * kappa.hat * I1k2) / (8 * pi * I0^2)
  R.fVM3 <- (4 * kappa.hat * I1k2 + 30 * kappa.hat^2 * I2k2 + 15 * kappa.hat^3 * I3k2) / (16 * pi * I0^2)
  R.fVM4 <- (8 * kappa.hat^2 * I0k2 + 105 * kappa.hat^4 * I2k2 +
               105 * kappa.hat^3 * I3k2 + 244 * kappa.hat^2 * I2k2) / (32 * pi * I0^2)

  R.hat <- 25 * R.fVM2 / 144 - 5 * R.fVM3 / 36 + R.fVM4 / 36
  bw <- (16 * sqrt(pi) / 3 * R.hat * n)^(2/9)

  return(bw)
}
