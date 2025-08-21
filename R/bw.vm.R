#' Compute the Optimal Bandwidth for Circular Data using von Mises Method
#'
#' This function computes the optimal smoothing parameter (bandwidth) for circular data
#' using the von Mises method. The optimal bandwidth is derived based on the maximum 
#' likelihood estimate of the concentration parameter kappa from a von Mises distribution
#' fitted to the data.
#'
#' @param x Data from which the smoothing parameter is to be computed. The object is
#'   coerced to a numeric vector in radians using `circular::conversion.circular`.
#'   Can be a numeric vector or an object of class `circular`.
#'
#' @return The computed optimal smoothing parameter, a numeric value derived from 
#'   the von Mises distribution fit to the data.
#'
#' @export
#'
#' @examples
#' # Example with circular data
#' library(circular)
#' set.seed(123)
#' x <- rvonmises(100, mu = circular(0), kappa = 2)
#' bw <- bw.vm(x)
#' print(bw)
#'
#' x <- rwrappednormal(100, mu = circular(1), rho = 0.7)
#' bw <- bw.vm(x)
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
bw.vm <- function(x) {
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
    cli::cli_alert_warning("{.var x} contains missing values, which will be removed.")
    x <- x[!is.na(x)]
  }
  n <- length(x)
  kappa.hat <- circular::mle.vonmises(x)$kappa
  
  b0.kappa <- besselI(kappa.hat, 0)
  b1.kappa <- besselI(kappa.hat, 1)
  b2.2kappa <- besselI(2 * kappa.hat, 2)
  
  # R_hat(f_VM^(2))
  r.fVM2 <- (3 * kappa.hat^2 * b2.2kappa + 2 * kappa.hat * b1.kappa) / (8 * pi * b0.kappa^2)

  kappa <- (2 * sqrt(pi) * r.fVM2 * n)^(2/5)
  return(kappa)
}
