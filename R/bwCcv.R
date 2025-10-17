#' @title Complete Cross-Validation for Circular Data
#'
#' @description This function calculates the optimal smoothing parameter (bandwidth) for circular data
#' using the complete cross-validation (CCV) method (see \doi{10.59170/stattrans-2024-024}). 
#'
#' @param x Data from which the smoothing parameter is to be computed. The object is
#'   coerced to a numeric vector in radians using \code{\link[circular]{circular}}.
#'   Can be a numeric vector or an object of class \code{circular}.
#' @param lower Lower boundary of the interval to be used in the search for the
#'   smoothing parameter \code{kappa}. Must be a positive numeric value less than \code{upper}.
#'   Default is 0.
#' @param upper Upper boundary of the interval to be used in the search for the
#'   smoothing parameter \code{kappa}. Must be a positive numeric value greater than \code{lower}.
#'   Default is 60.
#' @param tol Convergence tolerance for the \code{\link[stats]{optimize}} function, determining the
#'   precision of the optimization process. Also used to detect convergence near boundaries:
#'   if the optimal value is within \code{tol} of \code{lower} or \code{upper}, a warning
#'   is issued suggesting interval adjustment. Default is 0.1.
#' 
#' @details The complete cross-validation (CCV) method is an alternative for bandwidth 
#' selection, originally proposed by Jones (1991) for linear densities. Its adaptation 
#' to the circular setting was recently studied by Hasilová et al. (2024).
#' 
#' The method uses functionals \eqn{T_m} defined as:
#' \deqn{T_m(\kappa) = \frac{(-1)^m}{n(n-1)}\sum_{i=1}^n\sum_{j \neq i}^n K_{\kappa}^{(2m)}(\theta_{i} - \theta_{j})}
#' where \eqn{K_{\kappa}^{(2m)}} is the (2m)-th derivative of \eqn{K_{\kappa}}.
#' 
#' The CCV criterion can be expressed as:
#' \deqn{CCV(\kappa) = R(f(\kappa)) - T_0(\kappa) + \frac{1}{2}\bar{\sigma}_h^2 T_1(\kappa) + \frac{1}{24}(\eta_{2}^4(K_{\kappa}) - \eta_{4}(K_{\kappa}))T_2(\kappa)}
#' where \eqn{\eta_{j}(K_{\kappa})} denotes the j-th moment of the kernel.
#' 
#' For the von Mises kernel, the explicit CCV criterion becomes:
#' \deqn{CCV(\kappa) = \frac{1}{n^2} \sum_{i=1}^n \sum_{j=1}^n (K_{\kappa} * K_{\kappa})(\theta_i - \theta_j) - T_0(\kappa) + \frac{A_1(\kappa)}{2\kappa}T_1(\kappa) + \frac{2A_1^2(\kappa) - A_2(\kappa)}{8\kappa^2}T_2(\kappa)}
#' where \eqn{A_k(\kappa) = I_k(\kappa)/I_0(\kappa)} is the ratio of modified Bessel functions.
#' 
#' The optimal bandwidth is obtained by minimizing this criterion over the interval 
#' \code{[lower, upper]}.
#'
#' @return The computed optimal smoothing parameter \code{kappa}, a numeric concentration 
#' parameter (analogous to inverse radians) that minimizes the smoothed cross-validation 
#' criterion within the interval \code{[lower, upper]} and the value of objective function 
#' at that point. Higher values indicate sharper, more concentrated kernels and less 
#' smoothing; lower values indicate broader kernels and more smoothing. If the 
#' optimization fails, a warning is issued. 
#'
#' @export
#'
#' @examples
#' # Example with circular data
#' library(circular)
#' set.seed(123)
#' x <- rwrappednormal(100, mu = circular(2), rho = 0.5)
#' bw <- bwCcv(x)
#' print(round(bw$minimum, 2))
#'
#' x <- rvonmises(100, mu = circular(0), kappa = 1)
#' bw <- bwCcv(x)
#' print(round(bw$minimum, 2))
#'
#' @references
#' Hasilová, K., Horová, I., Valis, D., & Zámečník, S. (2024).
#' A comprehensive exploration of complete cross-validation for circular data.
#' \emph{Statistics in Transition New Series}, 25(3):1--12. \doi{10.59170/stattrans-2024-024}
#'
#' @seealso \link{bwScv}, \link{bwLscv}, \link{bwCcv}
#'
#' @importFrom stats optimize
#' @importFrom circular circular
#' @import cli
bwCcv <- function(x,
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
  x <- as.numeric(x)
  if (any(is.na(x))) {
    cli::cli_alert_warning("{.var x} contains missing values, which will be removed.")
    x <- x[!is.na(x)]
  }
  if (!is.numeric(lower)) {
    cli::cli_alert_warning(
      c(
        "Argument {.var lower} must be numeric. ",
        "Default value 0 for lower boundary was used."
      )
    )
    lower <- 0
  }
  if (!is.numeric(upper)) {
    cli::cli_alert_warning(
      c(
        "Argument {.var upper} must be numeric. ",
        "Default value 60 for upper boundary was used."
      )
    )
    upper <- 60
  }
  if (lower < 0 | lower >= upper) {
    cli::cli_alert_warning(
      c(
        "The boundaries must be positive numbers and 'lower' must be smaller than 'upper'. ",
        "Default boundaries lower=0, upper=60 were used."
      )
    )
    lower <- 0
    upper <- 60
  }

  ccv <- function(x, kappa) {
    kappa_min <- sqrt(.Machine$double.eps)
    if (kappa < kappa_min) {
      kappa <- kappa_min
    }
    
    n <- length(x)
    grid <- outer(x, x, "-")

    cos_grid <- cos(grid)
    sin_grid <- sin(grid)
    exp_kappa_cos <- exp(kappa * cos_grid)

    cos_0 <- 1
    sin_0 <- 0
    exp_0 <- exp(kappa * cos_0)

    b0_kappa <- suppressWarnings(besselI(kappa, 0))
    b0_kappacosgrid <- suppressWarnings(besselI(kappa * sqrt(2 * (1 + cos_grid)), 0))
    b1_kappa <- suppressWarnings(besselI(kappa, 1))
    b2_kappa <- suppressWarnings(besselI(kappa, 2))

    factor_1 <- 1 / (2 * pi * n ^ 2 * b0_kappa ^ 2)
    part_1 <- factor_1 * sum(b0_kappacosgrid)

    factor_2 <- 1 / (n * (n - 1) * 2 * pi * b0_kappa)
    part_2 <- factor_2 * (sum(exp_kappa_cos) - n * exp_0)

    factor_3 <- (1 / (2 * kappa)) * (b1_kappa / b0_kappa) * (1 / (2 * pi * b0_kappa * n * (n - 1)))
    arg_3 <- exp_kappa_cos * (kappa ^ 2 * sin_grid ^ 2 - kappa * cos_grid)
    arg_3_1 <- exp_0 * (-kappa)
    part_3 <- -factor_3 * (sum(arg_3) - n * arg_3_1)

    factor_4 <- (1 / (8 * kappa ^ 2)) * (2 * (b1_kappa / b0_kappa) ^ 2 - b2_kappa / b0_kappa) * (1 / (2 * pi * b0_kappa * n * (n - 1)))
    arg_4 <- exp_kappa_cos * (
      kappa ^ 4 * sin_grid ^ 4 - 6 * kappa ^ 3 * sin_grid ^ 2 * cos_grid +
        3 * kappa ^ 2 * (cos_grid ^ 2 - sin_grid ^ 2) - kappa ^
        2 * sin_grid ^ 2 + kappa * cos_grid
    )
    arg_4_1 <- exp_0 * (3 * kappa ^ 2 + kappa)
    part_4 <- factor_4 * (sum(arg_4) - n * arg_4_1)

    result <- part_1 - part_2 + part_3 + part_4
    return(result)
  }

  bw <- optimize(
    ccv,
    interval = c(lower, upper),
    maximum = FALSE,
    tol = tol,
    x = x
  )
  if (bw$minimum < lower + tol | bw$minimum > upper - tol) {
    cli::cli_alert_warning("Minimum/maximum occurred at one end of the range.")
  }
  return(bw)
}