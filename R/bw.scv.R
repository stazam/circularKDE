#' Compute the Optimal Bandwidth for Circular Data using Smoothed Cross-Validation
#'
#' This function computes the optimal smoothing parameter (bandwidth) for circular data
#' using a smoothed cross-validation (SCV) method (see <doi:doi.org/10.1007/s00180-023-01401-0>).
#' It searches for the value of the smoothing parameter `nu` that minimizes the SCV criterion within the
#' specified interval `[lower, upper]`.
#'
#' @param x Data from which the smoothing parameter is to be computed. The object is
#'   coerced to a numeric vector in radians using `circular::conversion.circular`.
#'   Can be a numeric vector or an object of class `circular`.
#' @param np An integer specifying the number of points used in numerical integration
#'   to evaluate the SCV criterion. A higher number increases precision but also
#'   computational cost. Default is 75.
#' @param lower Lower boundary of the interval for the optimization of the smoothing
#'   parameter `mu`. Must be a positive numeric value smaller than `upper`.
#'   Default is 0.
#' @param upper Upper boundary of the interval for the optimization of the smoothing
#'   parameter `mu`. Must be a positive numeric value greater than `lower`.
#'   Default is 60.
#' @param tol Convergence tolerance used in the `optimize` function. Determines how
#'   precisely the optimal value is estimated. Default is 0.1.
#'
#' @return The computed optimal smoothing parameter `mu`, the numeric value
#' that minimizes the smoothed cross-validation criterion.
#'
#' @export
#'
#' @examples
#' # Example with numeric data in radians
#' set.seed(123)
#' x <- runif(100, 0, 2 * pi)
#' bw <- bw.scv(x)
#' print(bw)
#'
#' # Example with circular data
#' library(circular)
#' x_circ <- rvonmises(100, mu = circular(0), kappa = 2)
#' bw <- bw.scv(x_circ)
#' print(bw)
#'
#' @importFrom circular conversion.circular
#' @importFrom stats optimize
#' @import cli

bw.scv <- function(x,
                   np = 75,
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
      c("{.var x} must be a numeric vector", "x" = "You've supplied a {.cls {class(x)}} vector.")
    )
  }
  x <- conversion.circular(
    x,
    units = "radians",
    zero = 0,
    rotation = "counter",
    modulo = "2pi"
  )
  attr(x, "class") <- attr(x, "circularp") <- NULL
  if (any(is.na(x))) {
    cli::cli_alert_warning("{.var x} contains missing values, which will be removed")
    x <- x[!is.na(x)]
  }
  if (!is.numeric(np)) {
    cli::cli_alert_warning(
      c(
        "Argument {.var np} must be numeric.",
        "Default value 75 for number of points for evalutaion of numerical integration was used."
      )
    )
    np <- 75
  }
  if (!is.numeric(lower)) {
    cli::cli_alert_warning(
      c(
        "Argument {.var lower} must be numeric.",
        "Default value 0 for lower boundary was used."
      )
    )
    upper <- 0
  }
  if (!is.numeric(upper)) {
    cli::cli_alert_warning(
      c(
        "Argument {.var upper} must be numeric.",
        "Default value 60 for upper boundary was used."
      )
    )
    upper <- 60
  }
  if (lower < 0 | lower >= upper) {
    cli::cli_alert_warning(
      c(
        "The boundaries must be positive numbers and 'lower' must be smaller that 'upper'",
        "Default boundaries lower=0, upper=60 were used"
      )
    )
    lower <- 0
    upper <- 60
  }
  bw <- optimize(
    scv,
    interval = c(lower, upper),
    maximum = FALSE,
    np = np,
    x = x
  )$minimum
  if (bw < lower + tol | bw > upper - tol) {
    cli::cli_alert_warning("Minimum/maximum occurred at one end of the range")
  }
  return(bw)
}



#              FISRT PART
trapezoidal.rule <- function(f,m,thetas,mu){

  h <- 2*pi / m
  knots <- seq(0,2 * pi, length = m)

  val <- 2*f(mu,knots,thetas)
  val <- sum(val) - val[1] - val[m]

  return ( 1/2 * h * val )
}


first <- function(mu,theta,thetas){

  first <- besselI(mu*sqrt(2*(1+ cos(theta - thetas[1]))),0)
  second <- besselI(mu*sqrt(2*(1+ cos(theta - thetas[2]))),0)

  return ( first * second )
}



first.part <- function(thetas, m, mu){

  n <- length(thetas)
  factor <- 1 /  (4 * pi^2 * n * (n - 1) * besselI(mu,0)^4)
  mat <- expand.grid(thetas,thetas)
  mat <-  apply(mat, 1,trapezoidal.rule, f = first, m = m, mu  = mu )
  mat[seq(1,length(mat),by = n+1)] <- 0
  return(factor * sum(mat))
}


#              SECOND PART
second <- function(mu,theta,thetas){


  first_part <- besselI(mu*sqrt(2*(1+ cos(theta - thetas[1]))),0)
  second_part <- exp(mu*cos(theta - thetas[2]))

  return ( first_part * second_part )

}



second.part <- function(thetas,m,mu){

  n <- length(thetas)

  factor <- 1/(2*n*(n-1)*pi^2 * besselI(mu,0)^3)

  mat <- expand.grid(thetas,thetas)
  mat <-  apply(mat, 1,trapezoidal.rule, f = second, m = m, mu  = mu )
  mat[seq(1,length(mat),by = n+1)] <- 0
  return(  factor *  sum(mat)  )

}


#              THIRD PART
third.part <- function(thetas,mu){

  n <- length(thetas)

  factor <- 1 / ( 2 * n * (n - 1) * pi * besselI(mu,0)^2)
  part <- besselI(mu * sqrt(2 * (1 + cos(outer(thetas, thetas, '-')))),0)

  return(  factor *  sum(part- diag(diag(part)))  )

}


scv <- function(mu,np,x){

  n <- length(x)
  first.part <- x %>% first.part(np,mu)
  second.part <- x %>% second.part(np,mu)
  third.part <- x %>% third.part(mu)

  return((first.part - second.part + third.part) + besselI(2*mu,0)/(2*n*pi*besselI(mu,0)^2))
}

