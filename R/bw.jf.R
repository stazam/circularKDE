#' @title Plug-in Method by Tsuruta and Sagae with additive Jones-Foster approach
#'
#' @description This function computes the optimal smoothing parameter (bandwidth) for circular data
#' using the plug-in method introduced by Tsuruta and Sagae (see \doi{10.1016/j.spl.2017.08.003}) with the
#' additive method from Jones and Foster (1993) to form higher-order kernel functions.
#'
#' @param x Data from which the smoothing parameter is to be computed. The object is
#'   coerced to a numeric vector in radians using \code{\link[circular]{conversion.circular}}.
#'   Can be a numeric vector or an object of class \code{circular}.
#' 
#' @details The plug-in approach estimates the optimal bandwidth through the following steps:
#' \enumerate{
#'   \item Apply the additive method from Jones and Foster (1993) to construct a p-th order kernel function.
#'   \item Plug-in the required functionals into AMISE expression.
#'   \item Solving for the bandwidth that minimizes the asymptotic mean integrated squared error (AMISE). The optimal
#'         bandwidth for the additive Jones-Foster method is given by:
#'         \deqn{\hat{\kappa}_{JF} = \left[\frac{16\sqrt{\pi}}{3} \hat{R}_{\hat{\tau}}\left(\frac{5f_{VM}^{(2)} + 2f_{VM}^{(4)}}{12}\right)n\right]^{2/9}}
#'         where the functional \eqn{\hat{R}_{\hat{\tau}}} is computed as a weighted linear combination under the von Mises assumption
#'         and \eqn{\hat{\tau}} is the MLE estimate of the von Mises concentration parameter used as the initial value.
#' }
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
#' \doi{10.1016/j.spl.2017.08.003}
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
    modulo = "2pi",
    template = "none"
  )
  attr(x, "class") <- attr(x, "circularp") <- NULL
  if (any(is.na(x))) {
    cli::cli_alert_warning("{.var x} contains missing values, which will be removed.")
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
