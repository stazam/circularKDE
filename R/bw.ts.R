#' Compute the Optimal Bandwidth for Circular Data using Taylor Series Method
#'
#' This function computes the optimal smoothing parameter (bandwidth) for circular data
#' using the Taylor Series method. The method is based on Taylor series approximations
#' to the von Mises functional and provides an alternative approach for bandwidth
#' selection in circular kernel density estimation.
#'
#' @param x Data from which the smoothing parameter is to be computed. The object is
#'   coerced to a numeric vector in radians using `circular::conversion.circular`.
#'   Can be a numeric vector or an object of class `circular`.
#'
#' @return The computed optimal smoothing parameter, a numeric value derived from
#'   the Taylor Series method for circular kernel density estimation.
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
#' García-Portugués, E. (2013). Exact risk improvement of bandwidth selectors
#' for kernel density estimation with directional data. \emph{Electronic
#' Journal of Statistics}, 7:1655--1685.
#' \doi{10.1214/13-ejs821}
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

  n <- length(x)
  nu.hat <- circular::mle.vonmises(x)$kappa
  b0.nu <- besselI(nu.hat, 0)
  b0.2nu <- besselI(2 * nu.hat, 0)
  b1.2nu <- besselI(2 * nu.hat, 1)
  b2.2nu <- besselI(2 * nu.hat, 2)
  b3.2nu <- besselI(2 * nu.hat, 3)

  r.fVM2 <- (3 * nu.hat^2 * b2.2nu + 2 * nu.hat * b1.2nu) / (8 * pi * b0.nu^2)
  r.fVM3 <- (4 * nu.hat * b1.2nu + 30 * nu.hat^2 * b2.2nu + 15 * nu.hat^3 * b3.2nu) / (16 * pi * b0.nu^2)
  r.fVM4 <- (8 * nu.hat^2 * b0.2nu + 105 * nu.hat^4 * b2.2nu +
               105 * nu.hat^3 * b3.2nu + 244 * nu.hat^2 * b2.2nu) / (32 * pi * b0.nu^2)

  # Approximation to m_vm functional (Appendix E, simplified form)
  R.hat <- 0.25 * r.fVM2 - 1.25 * r.fVM3 + 0.25 * r.fVM4

  bw <- (288 / (33 - 16 * sqrt(2) / sqrt(5)) * R.hat * n)^(2/9)

  return(bw)
}
