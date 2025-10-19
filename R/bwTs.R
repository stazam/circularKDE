#' @title Plug-in Method by Tsuruta and Sagae with multiplicative method from Terrell and Scott
#'
#' @description This function computes the optimal smoothing parameter (bandwidth) for circular data
#' using the plug-in method, introduced by Tsuruta and Sagae (see \doi{10.1016/j.spl.2017.08.003}) with the
#' multiplicative method from Terrell and Scott (1980) to form higher-order kernel functions.
#' The method optimally balances the bias-variance tradeoff by minimizing asymptotic mean
#' integrated squared error.
#'
#' @param x Data from which the smoothing parameter is to be computed. The object is
#'   coerced to a numeric vector in radians using \code{\link[circular]{conversion.circular}}.
#'   Can be a numeric vector or an object of class \code{circular}.
#' @param verbose Logical indicating whether to print intermediate computational values
#'   for debugging purposes. Useful for diagnosing floating-point instability in
#'   Bessel-based moment terms. Default is FALSE.
#'
#' @details The plug-in approach estimates the optimal bandwidth through the following steps:
#' \enumerate{
#'   \item Apply the multiplicative method from Terrell and Scott (1980) to construct a p-th order kernel function.
#'   \item Derive expression for asymptotic mean integrated squared error (AMISE) expression.
#'   \item Solving for the bandwidth that minimizes the AMISE. The optimal
#'         bandwidth for the multiplicative Terrell-Scott method is given by:
#'         \deqn{\hat{\kappa}_{TS} = \left[\frac{288}{33 - 16\sqrt{2/5}} \hat{R}_{\hat{\tau}}(m_{VM}) n\right]^{2/9}}
#'         where the computational formula is:
#'         \deqn{m_{VM}(\theta) := [2\{f_{VM}^{(2)}\}^2/f_{VM} - 5f_{VM}^{(2)} + 2f_{VM}^{(4)}]/4}
#'         and \eqn{\hat{R}_{\hat{\tau}}(m_{VM})} is the functional computed under the von Mises assumption
#'         using the multiplicative approach. The parameter \eqn{\hat{\tau}} is the MLE estimate of the
#'         von Mises concentration parameter used as the initial value.
#' }
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
#' bw <- bwTs(x)
#' print(bw)
#'
#' x <- rwrappednormal(100, mu = circular(1), rho = 0.7)
#' x <- as.numeric(x)
#' bw <- bwTs(x)
#' print(bw)
#'
#' @references
#' Tsuruta, Yasuhito & Sagae, Masahiko (2017). Higher order kernel density
#' estimation on the circle. \emph{Statistics & Probability Letters}, 131:46--50.
#' \doi{10.1016/j.spl.2017.08.003}
#'
#' Terrell, George R. & Scott, David W. (1980). On improving convergence rates
#' for nonnegative kernel density estimators. \emph{The Annals of Statistics},
#' 8(5):1160--1163.
#'
#' @seealso \link{bwScv}, \link{bwLscvg}, \link{bwCcv}
#'
#' @importFrom stats optimize
#' @importFrom circular mle.vonmises circular
#' @import cli
bwTs <- function(x, verbose = FALSE) {
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
    if (length(x) < 5) {
      cli::cli_abort(
        c("{.var x} must be a numeric vector of length at least 5 after removing missing values.",
          "x" = "You've supplied an object of length {length(x)} after removing missing values.")
      )
    }
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
  b2_kappa <- besselI(kappa_hat, 2)
  b2_2kappa <- besselI(2 * kappa_hat, 2)
  b3_2kappa <- besselI(2 * kappa_hat, 3)

  r_fVM2 <- (3 * kappa_hat^2 * b2_2kappa + 2 * kappa_hat * b1_2kappa) / (8 * pi * b0_kappa^2)
  r_fVM3 <- (4 * kappa_hat * b1_2kappa + 30 * kappa_hat^2 * b2_2kappa + 15 * kappa_hat^3 * b3_2kappa) / (16 * pi * b0_kappa^2)
  r_fVM4 <- (8 * kappa_hat^2 * b0_2kappa + 105 * kappa_hat^4 * b2_2kappa +
               105 * kappa_hat^3 * b3_2kappa + 244 * kappa_hat^2 * b2_2kappa) / (32 * pi * b0_kappa^2)
  r_fVM2VM1 <- (41 * kappa_hat^4 * b2_2kappa + 12 * kappa_hat^2 * b2_kappa - 87 * kappa_hat^3 * b3_2kappa) / (32 * pi * b0_kappa^2)
  integral_1 <- (4 * kappa_hat^3 * b1_2kappa - 14 * kappa_hat^2 * b2_2kappa - 3 * kappa_hat^3 * b3_2kappa) / (16 * pi * b0_kappa^2)
  integral_2 <- (36 * kappa_hat^2 * b2_2kappa + 9 * kappa_hat^4 * b2_2kappa + 25 * kappa_hat^3 * b3_2kappa) / (32 * pi * b0_kappa^2)

  if (verbose) {
    cli::cli_alert_info("Moment terms:")
    cli::cli_alert_info("  r_fVM2 = {.val {round(r_fVM2, 8)}}")
    cli::cli_alert_info("  r_fVM3 = {.val {round(r_fVM3, 8)}}")
    cli::cli_alert_info("  r_fVM4 = {.val {round(r_fVM4, 8)}}")
    cli::cli_alert_info("  r_fVM2VM1 = {.val {round(r_fVM2VM1, 8)}}")
    cli::cli_alert_info("  integral_1 = {.val {round(integral_1, 8)}}")
    cli::cli_alert_info("  integral_2 = {.val {round(integral_2, 8)}}")
  }

  r_hat <- 1/4 * r_fVM2VM1 + 25/16 * r_fVM2 - 5/4 * r_fVM3 + 1/4 * r_fVM4 - 5/4 * integral_1 - 1/2 * integral_2

  if (!is.finite(r_hat) || r_hat <= 0) {
    cli::cli_abort(
      c("Computed functional R-hat is invalid for bandwidth calculation.",
        "x" = "Computed r_hat: {.val {r_hat}}",
        "i" = "This may indicate numerical instability or inappropriate data distribution.")
    )
  }

  c_ts <- 288 / (33 - 16 * sqrt(2/5))
  bw <- (c_ts * r_hat * n)^(2/9)

  return(bw)
}



