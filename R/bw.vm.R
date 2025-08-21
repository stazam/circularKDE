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
#' García-Portugués, E. (2013). Exact risk improvement of bandwidth selectors
#' for kernel density estimation with directional data. \emph{Electronic
#' Journal of Statistics}, 7:1655--1685.
#' \doi{10.1214/13-ejs821}
#'
#' @importFrom stats optimize
#' @import circular
#' @import cli
bw.vm <- function(x) {
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
  n <- length(x)
  kappa_hat <- circular::mle.vonmises(x)$kappa
  
  I0 <- besselI(kappa_hat, 0)
  I1 <- besselI(kappa_hat, 1)
  I2 <- besselI(2 * kappa_hat, 2)
  
  # R_hat(f_VM^(2))
  R_fVM2 <- (3 * kappa_hat^2 * I2 + 2 * kappa_hat * I1) / (8 * pi * I0^2)
  
  kappa_VM <- (2 * sqrt(pi) * R_fVM2 * n)^(2/5)
  return(kappa_VM)
}
