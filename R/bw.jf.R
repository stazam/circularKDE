#' Compute the Optimal Bandwidth for Circular Data using circular version of the additive method from Jones and Foster.
#'
#' This function computes the optimal smoothing parameter (bandwidth) for circular data
#' using the circular version of the additive method from Jones and Foster (1993). The method
#' is based on the new class of higher order kernel functions with new moments introduced by Tsurunga (see <doi:doi.org/10.1016/j.spl.2017.08.003>).
#'
#' @param x Data from which the smoothing parameter is to be computed. The object is
#'   coerced to a numeric vector in radians using `circular::conversion.circular`.
#'   Can be a numeric vector or an object of class `circular`.
#'
#' @return The computed optimal smoothing parameter, a numeric value derived from
#'   the circular version of the additive method for circular kernel density estimation.
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
#' Tsuruta, Yasuhito & Sagae, Masahiko (2017). Higher order kernel density
#' estimation on the circle. \emph{Statistics & Probability Letters}, 131:46--50.
#' \doi{10.1016/j.spl.2017.07.027}
#'
#' @importFrom stats optimize
#' @import circular
#' @import cli
bw.jf <- function(x) {
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
  b2.2kappa <- besselI(2 * kappa.hat, 2)
  b3.2kappa <- besselI(2 * kappa.hat, 3)

  r.fVM2 <- (3 * kappa.hat^2 * b2.2kappa + 2 * kappa.hat * b1.2kappa) / (8 * pi * b0.kappa^2)
  r.fVM3 <- (4 * kappa.hat * b1.2kappa + 30 * kappa.hat^2 * b2.2kappa + 15 * kappa.hat^3 * b3.2kappa) / (16 * pi * b0.kappa^2)
  r.fVM4 <- (8 * kappa.hat^2 * b0.2kappa + 105 * kappa.hat^4 * b2.2kappa +
               105 * kappa.hat^3 * b3.2kappa + 244 * kappa.hat^2 * b2.2kappa) / (32 * pi * b0.kappa^2)

  r.hat <- 25 * r.fVM2 / 144 - 5 * r.fVM3 / 36 + r.fVM4 / 36
  bw <- ((16 * sqrt(pi)) / 3 * (r.hat * n))^(2/9)

  return(bw)
}
