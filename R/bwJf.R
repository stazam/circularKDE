#' @title Plug-in Method by Tsuruta and Sagae with additive Jones-Foster approach
#'
#' @description This function computes the optimal smoothing parameter (bandwidth) for circular data
#' using the plug-in method introduced by Tsuruta and Sagae (see \doi{10.1016/j.spl.2017.08.003}) with the
#' additive method from Jones and Foster (1993) to form higher-order kernel functions.
#'
#' @param x Data from which the smoothing parameter is to be computed. The object is
#'   coerced to a numeric vector in radians using \code{\link[circular]{circular}}.
#'   Can be a numeric vector or an object of class \code{circular}.
#' @param verbose Logical indicating whether to print intermediate computational values
#'   for debugging and teaching purposes. Shows kappa_hat, r_hat, and component
#'   calculations. Default is FALSE.
#'
#' @details The plug-in approach estimates the optimal bandwidth through the following steps:
#' \enumerate{
#'   \item Apply the additive method from Jones and Foster (1993) to construct a p-th order kernel function.
#'   \item Derive expression for asymptotic mean integrated squared error (AMISE) expression.
#'   \item Solving for the bandwidth that minimizes the AMISE. The optimal
#'         bandwidth for the additive Jones-Foster method is given by:
#'         \deqn{\hat{\kappa}_{JF} = \left[\frac{16\sqrt{\pi}}{3} \hat{R}_{\hat{\tau}}\left(\frac{5f_{VM}^{(2)} + 2f_{VM}^{(4)}}{12}\right)n\right]^{2/9}}
#'         where the functional \eqn{\hat{R}_{\hat{\tau}}} is computed as a weighted linear combination under the von Mises assumption
#'         and \eqn{\hat{\tau}} is the MLE estimate of the von Mises concentration parameter used as the initial value.
#' }
#'
#' @return The computed optimal smoothing parameter \code{kappa}, a numeric concentration
#'  parameter (analogous to inverse radians) derived from the circular version of the
#'  additive method for circular kernel density estimation. Higher values indicate sharper,
#'  more concentrated kernels and less smoothing; lower values indicate broader kernels
#'  and more smoothing.
#'
#' @export
#'
#' @examples
#' # Example with circular data
#' library(circular)
#' set.seed(123)
#' x <- rvonmises(100, mu = circular(0), kappa = 2)
#' bw <- bwJf(x)
#' print(bw)
#'
#' x <- rwrappednormal(100, mu = circular(1), rho = 0.7)
#' bw <- bwJf(x)
#' print(bw)
#'
#' @references
#' Tsuruta, Yasuhito & Sagae, Masahiko (2017). Higher order kernel density
#' estimation on the circle. \emph{Statistics & Probability Letters}, 131:46--50.
#' \doi{10.1016/j.spl.2017.08.003}
#'
#' Jones, M. C. & Foster, P. J. (1993). Generalized jackknifing and higher-order kernels.
#' \emph{Journal of Nonparametric Statistics}, 3:81--94.
#' \doi{10.1080/10485259308832573}
#'
#' @seealso \link{bwScv}, \link{bwLscvg}, \link{bwCcv}
#'
#' @importFrom stats optimize
#' @importFrom circular mle.vonmises circular
#' @import cli
bwJf <- function(x, verbose = FALSE) {
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
  if (any(is.na(x))) {
    cli::cli_alert_warning("{.var x} contains missing values, which will be removed.")
    x <- x[!is.na(x)]
  }

  n <- length(x)
  kappa_hat <- mle.vonmises(x)$kappa

  if (!is.finite(kappa_hat) || kappa_hat <= 0) {
    cli::cli_abort(
      c("MLE estimation of concentration parameter failed.",
        "x" = "Computed kappa_hat: {.val {kappa_hat}}",
        "i" = "This may indicate insufficient data or degenerate distribution.")
    )
  }


  b0_kappa <- besselI(kappa_hat, 0)
  b0_2kappa <- besselI(2 * kappa_hat, 0)
  b1_2kappa <- besselI(2 * kappa_hat, 1)
  b2_2kappa <- besselI(2 * kappa_hat, 2)
  b3_2kappa <- besselI(2 * kappa_hat, 3)

  r_fVM2 <- (3 * kappa_hat^2 * b2_2kappa + 2 * kappa_hat * b1_2kappa) / (8 * pi * b0_kappa^2)
  r_fVM3 <- (4 * kappa_hat * b1_2kappa + 30 * kappa_hat^2 * b2_2kappa + 15 * kappa_hat^3 * b3_2kappa) / (16 * pi * b0_kappa^2)
  r_fVM4 <- (8 * kappa_hat^2 * b0_2kappa + 105 * kappa_hat^4 * b2_2kappa +
               105 * kappa_hat^3 * b3_2kappa + 244 * kappa_hat^2 * b2_2kappa) / (32 * pi * b0_kappa^2)

  # formula based on approximation of R(f_{VM}) see suplementary material of Tsuruta and Sagae (2017)
  r_hat <- 25 * r_fVM2 / 144 - 5 * r_fVM3 / 36 + r_fVM4 / 36
  if (verbose) {
    cli::cli_alert_info("MLE concentration parameter estimate:")
    cli::cli_alert_info("  kappa_hat = {.val {round(kappa_hat, 6)}}")
    cli::cli_alert_info("Computed components:")
    cli::cli_alert_info("  r_hat = {.val {round(r_hat, 6)}}")
  }
  bw <- ((16 * sqrt(pi)) / 3 * (r_hat * n))^(2/9)

  return(bw)
}
