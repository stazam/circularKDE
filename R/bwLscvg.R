#' @title Generalized Least Squares Cross-Validation for Circular Data
#'
#' @description This function computes the optimal smoothing parameter (bandwidth) for circular data
#' using a generalized least squares cross-validation method (see \doi{10.1007/s00180-023-01401-0}).
#'
#' @param x Data from which the smoothing parameter is to be computed. The object is
#'   coerced to a numeric vector in radians using \code{\link[circular]{conversion.circular}}.
#'   Can be a numeric vector or an object of class \code{circular}.
#' @param g A numeric scalar that controls the variability in the cross-validation
#'   procedure. It influences the scaling in the internal calculations, affecting the
#'   bandwidth estimation. Default is 4.
#' @param lower Lower boundary of the interval to be used in the search for the
#'   smoothing parameter \code{kappa}. Must be a positive numeric value less than \code{upper}.
#'   Default is 0.
#' @param upper Upper boundary of the interval to be used in the search for the
#'   smoothing parameter \code{kappa}. Must be a positive numeric value greater than \code{lower}.
#'   Default is 60.
#' @param tol Convergence tolerance for the \code{\link[stats]{optimize}} function, determining the
#'   precision of the optimization process. Default is 0.1.
#' 
#' @details The generalized least squares cross-validation method (LSCV_g) is an 
#' adaptation of the method originally introduced by Zhang for linear data, developed 
#' for circular data (see Zamecnik, et.al 2025) to address the finite-sample performance issues of standard LSCV.
#' 
#' The LSCV_g criterion is defined as:
#' \deqn{LSCV_g(\kappa) = \frac{1}{n}R(K_{\kappa}) + \frac{1}{n(n-1)} \sum_{i=1}^n \sum_{j \neq i}^n \left(\frac{n-1}{n} (K_{\kappa}*K_{\kappa})(\theta_i-\theta_j) + \frac{2}{g(g-2)} K_{\kappa/g}(\theta_i-\theta_j) - \frac{g-1}{g-2} K_{\kappa/2}(\theta_i-\theta_j)\right)}
#' 
#' Using the von Mises kernel, this takes the computational form:
#' \deqn{LSCV_g(\kappa) = \frac{1}{2\pi n^2} \sum_{i=1}^n \sum_{j=1}^n \frac{I_0(\kappa \sqrt{2(1+\cos(\theta_i-\theta_j))})}{I_0^2(\kappa)} + \frac{1}{n(n-1)} \sum_{i=1}^n \sum_{j \neq i}^n \left(\frac{2}{g(g-2)} \frac{\exp(\frac{\kappa}{g}\cos(\theta_i-\theta_j))}{2\pi I_0(\kappa/g)} - \frac{g-1}{g-2} \frac{\exp(\frac{\kappa}{2}\cos(\theta_i-\theta_j))}{2\pi I_0(\kappa/2)}\right)}
#' 
#' The optimal bandwidth is obtained by minimizing this criterion over the interval 
#' \code{[lower, upper]}. 
#'
#' @return The computed optimal smoothing parameter \code{kappa}, a numeric value that
#'   minimizes the least squares cross-validation criterion within the interval
#'   \code{[lower, upper]}.
#'
#' @export
#'
#' @examples
#' # Example with circular data
#' library(circular)
#' set.seed(123)
#' x <- rwrappednormal(100, mu = circular(2), rho = 0.5)
#' bw <- bwLscvg(x)
#' print(bw)
#'
#' x <- rvonmises(100, mu = circular(0), kappa = 1)
#' bw <- bwLscvg(x)
#' print(bw)
#' 
#' @references
#' Zámečník, S., Horová, I., & Hasilová, K. (2025). Generalised least square 
#' cross-validation for circular data. \emph{Communications in Statistics}. 
#' \doi{10.1007/s00180-023-01401-0}
#' 
#' Zhang, J. (2015). Generalized least squares cross-validation in kernel density 
#' estimation. \emph{Statistica Neerlandica}, 69(3), 315-328. 
#' \doi{10.1111/stan.12061}
#'
#' @importFrom stats optimize
#' @import circular
#' @import cli
bwLscvg <- function(x,
                     g = 4,
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
    modulo = "2pi",
    template = "none"
  )
  x <- as.numeric(x)
  if (any(is.na(x))) {
    cli::cli_alert_warning("{.var x} contains missing values, which will be removed.")
    x <- x[!is.na(x)]
  }
  if (!is.numeric(g)) {
    cli::cli_alert_warning(c(
      "Argument {.var g} must be numeric. ",
      "Default value 4 for coefficient was used."
    ))
    g <- 4
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

  lscvg <- function(x, kappa, g) {
    n <- length(x)
    grid <- outer(x, x, '-')
    cos_grid <- cos(grid)

    b0_kappa      <- besselI(kappa, 0)
    b0_2_kappa     <- besselI(2 * kappa, 0)
    b0_kappa_g    <- besselI(kappa / g, 0)
    b0_kappa_half <- besselI(kappa / 2, 0)

    part_1 <- b0_2_kappa / (2 * pi * n * (b0_kappa ^ 2))

    factor_2 <- (2 / (g * (g - 2))) * (1 / (2 * pi * b0_kappa_g))
    part_2 <- factor_2 * (sum(exp((kappa / g) * cos_grid)) - n * exp(kappa / g))

    arg_3 <- kappa * sqrt(2 * (1 + cos_grid))
    factor_3 <- ((n - 1) / n) * (1 / (2 * pi * (b0_kappa ^ 2)))
    part_3 <- factor_3 * (sum(besselI(arg_3, 0)) - n * b0_2_kappa)

    factor_4 <- ((g - 1) / (g - 2)) * (1 / (2 * pi * b0_kappa_half))
    part_4 <-  factor_4 * (sum(exp((kappa / 2) * cos_grid)) - n * exp(kappa / 2))

    result <- part_1 + (part_2 + part_3 - part_4) / (n * (n - 1))
    return(result)
  }
  bw <- optimize(
    lscvg,
    interval = c(lower, upper),
    maximum = FALSE,
    x = x ,
    g = g
  )$minimum
  if (bw < lower + tol | bw > upper - tol) {
    cli::cli_alert_warning("Minimum/maximum occurred at one end of the range.")
  }
  return(bw)
}
