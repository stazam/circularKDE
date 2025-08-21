#' Compute the Optimal Bandwidth for Circular Data using Fourier Method
#'
#' This function computes the optimal smoothing parameter (bandwidth) for circular data
#' using the Fourier method. The method estimates the bandwidth based on Fourier
#' coefficients and utilizes a data-driven approach to select the optimal truncation
#' point for the Fourier series expansion.
#'
#' @param x Data from which the smoothing parameter is to be computed. The object is
#'   coerced to a numeric vector in radians using `circular::conversion.circular`.
#'   Can be a numeric vector or an object of class `circular`.
#'
#' @return The computed optimal smoothing parameter, a numeric value derived from
#'   the Fourier method for circular kernel density estimation.
#'
#' @export
#'
#' @examples
#' # Example with circular data
#' library(circular)
#' set.seed(123)
#' x <- rvonmises(100, mu = circular(0), kappa = 2)
#' bw <- bw.fo(x)
#' print(bw)
#'
#' x <- rwrappednormal(100, mu = circular(1), rho = 0.7)
#' bw <- bw.fo(x)
#' print(bw)
#'
#' @references
#' Tenreiro, C. (2022). Kernel density estimation for circular data: a Fourier 
#' series-based plug-in approach for bandwidth selection. \emph{Journal of 
#' Nonparametric Statistics}, 34(2):377--406.
#' \doi{10.1080/10485252.2022.2057974}
#'
#' @importFrom stats optimize
#' @import circular
#' @import cli
bw.fo <- function(x) {
  
  n <- length(x)
  C1 <- 0.25
  C2 <- 25      
  gamma <- 0.5

  Ln <- trunc(C1 * n^(1/11)) + 1
  Un <- trunc(C2 * n^(1/11)) 
  
  theta2 <- numeric(Un)
  H <- numeric(Un)
  
  c1 <- (mean(cos(x)))^2 + (mean(sin(x)))^2
  theta2[1] <- c1 / pi 
  
  c_bar <- (n / (n - 1)) * (c1 - 1/n)
  H[1] <- 1/n - gamma * (1 + 1/n) * c_bar
  
  for (k in 2:Un) {
    c_k <- (mean(cos(k * x)))^2 + (mean(sin(k * x)))^2
    theta2[k] <- theta2[k - 1] + (k^4 * c_k) / pi
    c_bar <- c_bar + (n / (n - 1)) * (c_k - 1/n)
    H[k] <- k / n - gamma * (1 + 1/n) * c_bar
  }
  m <- (Ln - 1) + which.min(H[Ln:Un])
  bandwidth <- (4 * pi)^(-1/10) * (theta2[m] * n)^(-1/5)
  return(bandwidth)
}