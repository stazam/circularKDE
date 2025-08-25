#' Compute the Optimal Bandwidth for Circular Data using circular version of multiplicative method from Terrell and Scott.
#'
#' This function computes the optimal smoothing parameter (bandwidth) for circular data
#' using the circular version of the multiplicative method from Terrell and Scott (1980). The method
#' is based on the new class of higher order kernel functions with new moments introduced by Tsurunga (see <doi:doi.org/10.1016/j.spl.2017.08.003>).
#'
#' @param x Data from which the smoothing parameter is to be computed. The object is
#'   coerced to a numeric vector in radians using `circular::conversion.circular`.
#'   Can be a numeric vector or an object of class `circular`.
#'
#' @return The computed optimal smoothing parameter, a numeric value derived from
#'   the circular version of the multiplicative method for circular kernel density estimation.
#'
#' @export
#'
#' @examples
#' # Example with circular data
#' library(circular)
#' set.seed(123)
#' x <- rvonmises(100, mu = circular(0), kappa = 2)
#' bw <- bw.ts(x)
#' print(bw)
#'
#' x <- rwrappednormal(100, mu = circular(1), rho = 0.7)
#' bw <- bw.ts(x)
#' print(bw)
#'
#' @references
#' Tsuruta, Yasuhito & Sagae, Masahiko (2017). Higher order kernel density
#' estimation on the circle. \emph{Statistics & Probability Letters}, 131:46--50.
#' \doi{10.1016/j.spl.2017.07.027}
#'
#' @importFrom stats optimize
#' @import circular
#' @import cli
bw.ts <- function(x) {
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
  if (any(is.na(x))) {
    cli::cli_alert_warning("{.var x} contains missing values, which will be removed")
    x <- x[!is.na(x)]
  }

  n <- length(x)
  kappa.hat <- circular::mle.vonmises(x)$kappa
  b0.kappa <- besselI(kappa.hat, 0)
  b0.2kappa <- besselI(2 * kappa.hat, 0)
  b1.2kappa <- besselI(2 * kappa.hat, 1)
  b2.kappa <- besselI(kappa.hat, 2)
  b2.2kappa <- besselI(2 * kappa.hat, 2)
  b3.2kappa <- besselI(2 * kappa.hat, 3)

  r.fVM2 <- (3 * kappa.hat^2 * b2.2kappa + 2 * kappa.hat * b1.2kappa) / (8 * pi * b0.kappa^2)
  r.fVM3 <- (4 * kappa.hat * b1.2kappa + 30 * kappa.hat^2 * b2.2kappa + 15 * kappa.hat^3 * b3.2kappa) / (16 * pi * b0.kappa^2)
  r.fVM4 <- (8 * kappa.hat^2 * b0.2kappa + 105 * kappa.hat^4 * b2.2kappa +
               105 * kappa.hat^3 * b3.2kappa + 244 * kappa.hat^2 * b2.2kappa) / (32 * pi * b0.kappa^2)
  r.fVM2VM1 <- (41 * kappa.hat^4 * b2.2kappa + 12 * kappa.hat^2 * b2.kappa - 87 * kappa.hat^3 * b3.2kappa) / (32 * pi * b0.kappa^2)
  integral.1 <- (4 * kappa.hat^3 * b1.2kappa - 14 * kappa.hat^2 * b2.2kappa - 3 * kappa.hat^3 * b3.2kappa) / (16 * pi * b0.kappa^2)
  integral.2 <- (36 * kappa.hat^2 * b2.2kappa + 9 * kappa.hat^4 * b2.2kappa + 25 * kappa.hat^3 * b3.2kappa) / (32 * pi * b0.kappa^2)

  r.hat <- 0.25 * r.fVM2VM1 + 1.5625 * r.fVM2 - 1.25 * r.fVM3 + 0.25 * r.fVM4 - 1.25 * integral.1 - 0.5 * integral.2

  bw <- (288 / (33 - 16 * sqrt(2) / sqrt(5)) * r.hat * n)^(2/9)

  return(bw)
}
