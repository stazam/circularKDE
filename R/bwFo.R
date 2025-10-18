#' @title Optimal Bandwidth Selection via Fourier Plug-in Method
#'
#' @description This function computes the optimal smoothing parameter (bandwidth) for circular data
#' using the Fourier series-based direct plug-in approach based on delta sequence estimators (see \doi{10.1080/10485252.2022.2057974}).
#'
#' @param x Data from which the smoothing parameter is to be computed. The object is
#'   coerced to a numeric vector in radians using \code{\link[circular]{circular}}.
#'   Can be a numeric vector or an object of class \code{circular}.
#' 
#' @details The Fourier-based plug-in estimator computes the optimal bandwidth using the formula:
#' \deqn{\hat{h}_{FO} := (4\pi)^{-1/10} \hat{\theta}_{2,\hat{m}}^{-1/5} n^{-1/5}}
#' where \eqn{\hat{\theta}_{2,\hat{m}}} is the estimator of the second-order functional 
#' \eqn{\theta_2(f)} based on the selected number of Fourier coefficients \eqn{\hat{m}}.
#' 
#' Under the assumption of von Mises density, this formula becomes:
#' \deqn{\hat{h}_{VM} = (4\pi)^{-1/10} \left(\frac{3\hat{\kappa}^2 I_0(2\hat{\kappa}) - \hat{\kappa}I_1(2\hat{\kappa})}{8\pi I_0(\hat{\kappa})^2}\right)^{-1/5} n^{-1/5}}
#' where \eqn{I_0} and \eqn{I_1} are the modified Bessel functions of the first kind of orders 0 and 1, 
#' and \eqn{\hat{\kappa}} is the estimated concentration parameter of the von Mises distribution.
#' 
#' @return The computed optimal smoothing parameter \code{kappa}, a numeric concentration
#'  parameter (analogous to inverse radians) derived from the Fourier method for circular 
#'  kernel density estimation. Higher values indicate sharper,  more concentrated kernels 
#'  and less smoothing; lower values indicate broader kernels and more smoothing.
#'
#' @export
#'
#' @examples
#' # Example with circular data
#' library(circular)
#' set.seed(123)
#' x <- rvonmises(100, mu = circular(0), kappa = 2)
#' bw <- bwFo(x)
#' print(bw)
#'
#' x <- rwrappednormal(100, mu = circular(1), rho = 0.7)
#' bw <- bwFo(x)
#' y <- density(x, bw=bw) 
#' plot(y, main="KDE with Fourier Plug-in Bandwidth")
#'
#' @references
#' Tenreiro, C. (2022). Kernel density estimation for circular data: a Fourier 
#' series-based plug-in approach for bandwidth selection. \emph{Journal of 
#' Nonparametric Statistics}, 34(2):377--406.
#' \doi{10.1080/10485252.2022.2057974}
#' 
#' @seealso \link{bwScv}, \link{bwLscvg}, \link{bwCcv}
#'
#' @importFrom stats optimize
#' @importFrom circular circular
#' @import cli
bwFo <- function(x, C1 = 0.25, C2 = 25, gamma = 0.5) {
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
  if (gamma <= 0 || !is.numeric(gamma) || length(gamma) != 1 || gamma > 1) {
    cli::cli_abort("{.var gamma} must be a strictly positive numeric value from the interval (0, 1].")
  }
  for (nm in c("C1", "C2")) {
    val <- get(nm, inherits = TRUE)
    if (!is.numeric(val) || length(val) != 1 || val <= 0) {
      cli::cli_abort("{.var {nm}} must be a single strictly positive numeric value.")
      }
  }
  if (C1 > C2) {
    cli::cli_abort("{.var C1} must be less than {.var C2}.")
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
    if (length(x) < 5) {
      cli::cli_abort(
        c("{.var x} must be a numeric vector of length at least 5 after removing missing values.", 
          "x" = "You've supplied an object of length {length(x)} after removing missing values.")
      )
    }
  }
  
  n <- length(x)
  # Ln, Un are the determined lower and upper bounds for m computed based on n a strictly positive constants C1 and C2
  # Choice of n^(1/11) provides the best rate of convergence for relative errors of the estiamtor
  Ln <- trunc(C1 * n^(1/11)) + 1
  Un <- trunc(C2 * n^(1/11)) 
  
  # Function H(m) to be minimized to select optimal m (the formula for H(m) incorporates c_bar)
  c1 <- (mean(cos(x)))^2 + (mean(sin(x)))^2  
  c_bar <- (n / (n - 1)) * (c1 - 1/n)

  H <- numeric(Un)
  H[1] <- 1/n - gamma * (1 + 1/n) * c_bar

  # theta2(m) - The only unknown part of fourier estimator is quadratic functional  
  theta2 <- numeric(Un)
  theta2[1] <- c1 / pi 
  
  # Compute H(k) and theta2(k) for k = 2, ..., Un
  for (k in 2:Un) {
    c_k <- (mean(cos(k * x)))^2 + (mean(sin(k * x)))^2
    theta2[k] <- theta2[k - 1] + (k^4 * c_k) / pi
    c_bar <- c_bar + (n / (n - 1)) * (c_k - 1/n)
    H[k] <- k / n - gamma * (1 + 1/n) * c_bar
  }
  # Optimal number of Fourier terms m
  m <- (Ln - 1) + which.min(H[Ln:Un])
  # Fourier plug-in bandwidth
  bw <- (4 * pi)^(-1/10) * (theta2[m] * n)^(-1/5)
  return(bw)
}