#' Title
#'
#' @param x Data from which the smoothing parameter is to be computed. The object is coerced to class circular.
#' @param g parameter for bigger variability
#' @param lower Lower boundary of the interval to be used in the search for the value of the smoothing parameter. Default value lower=0.
#' @param upper Upper boundary of the interval to be used in the search for the value of the smoothing parameter. Default value upper=60.
#' @param tol Convergence tolerance for optimize.
#'
#' @return something
#' @export
#'
#' @examples something
#'
#' @importFrom circular conversion.circular
#' @importFrom stats optimize
#' @import cli
bw.lscvg <- function(x,
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
      c("{.var x} must be a numeric vector", "x" = "You've supplied a {.cls {class(x)}} vector.")
    )
  }
  x <- conversion.circular(
    x,
    units = "radians",
    zero = 0,
    rotation = "counter",
    modulo = "2pi",
    template = "none"
  )
  attr(x, "class") <- attr(x, "circularp") <- NULL
  if (any(is.na(x))) {
    cli::cli_alert_warning("{.var x} contains missing values, which will be removed")
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
        "The boundaries must be positive numbers and 'lower' must be smaller that 'upper'. ",
        "Default boundaries lower=0, upper=60 were used."
      )
    )
    lower <- 0
    upper <- 60
  }

  lscvg <- function(x, nu, g) {
    n <- length(x)
    grid <- outer(x, x, '-')

    part1 <- 1 / n * besselI(2 * nu, 0) / (2 * pi * besselI(nu, 0) ^ 2)

    factor2.1 <- 2 / (g * (g - 2))
    factor2.2 <-  1 / (2 * pi * besselI(nu / g, 0))
    part2 <- factor2.1 * factor2.2 * (sum(exp(nu / g * cos(grid))) - sum(exp(nu / g * cos(rep(
      0, n
    )))))

    factor3.1 <- (n - 1) / n
    factor3.2 <-  1 / (2 * pi * besselI(nu, 0) ^ 2)
    part3 <- factor3.1 * factor3.2 * (sum(besselI(nu * sqrt(2 * (
      1 + cos(grid)
    )), 0)) - sum(besselI(nu * sqrt(2 * (
      1 + cos(rep(0, n))
    )), 0)))


    factor4.1 <- (g - 1) / (g - 2)
    factor4.2 <-  1 / (2 * pi * besselI(nu / 2, 0))
    part4 <- factor4.1 * factor4.2 * (sum(exp(nu / 2 * cos(grid))) - sum(exp(nu / 2 * cos(rep(
      0, n
    )))))

    output <- part1 +  1 / (n * (n - 1)) * (part2 + part3 - part4)
    return(output)
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




