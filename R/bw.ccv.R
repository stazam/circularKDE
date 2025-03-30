#' Title
#'
#' @param x Data from which the smoothing parameter is to be computed. The object is coerced to class circular.
#' @param lower Lower boundary of the interval to be used in the search for the value of the smoothing parameter. Default value lower=0.
#' @param upper Upper boundary of the interval to be used in the search for the value of the smoothing parameter. Default value upper=60.
#' @param tol Convergence tolerance for optimize.
#'
#' @return something
#' @export
#'
#' @examples something
bw.ccv <- function(x,
                   lower = 0,
                   upper = 60,
                   tol = 0.1) {
  n <- length(x)
  if (n == 0) {
    cli::cli_abort(
      c("{.var x} must be a non-empty object.", "x" = "You've supplied an object of length {n}.")
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
  if (!is.numeric(x)) {
    cli::cli_abort(
      c("{.var x} must be a numeric vector", "x" = "You've supplied a {.cls {class(x)}} vector.")
    )
  }
  if (any(is.na(x))) {
    cli::cli_alert_warning("{.var x} contains missing values, which will be removed")
    x <- x[!is.na(x)]
  }
  n <- length(x)
  if (n == 0) {
    cli_abort("{.var x} is after removal of length {.var n}.", )
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

  ccv <- function(x, nu) {
    n <- length(x)
    grid <- outer(x, x, '-')

    factor2.1 <- 1  / (2 * pi * n ^ 2 * besselI(nu, 0) ^ 2)
    part2 <- factor2.1 * sum(besselI(nu * sqrt(2 * (1 + cos(
      grid
    ))), 0))


    factor3.1 <-  1 / (n * (n - 1))
    factor3.2 <- 1 / (2 * pi * besselI(nu, 0))
    part3 <- factor3.1 * factor3.2 *  (sum(exp(nu * cos(grid))) - sum(exp(nu * cos(rep(
      0, n
    )))))


    factor4.1 <- 1 / (2 * nu)
    factor4.2 <- besselI(nu, 1) / besselI(nu, 0)
    factor4.3 <- 1 / (2 * pi * besselI(nu, 0) * n * (n - 1))
    argument4.1 <- exp(nu * cos(grid)) * (nu ^ 2 * sin(grid) ^ 2 - nu * cos(grid))
    argument4.2 <- exp(nu * cos(rep(0, n))) * (nu ^ 2 * sin(rep(0, n)) ^
                                                 2 - nu * cos(rep(0, n)))
    part4 <- -factor4.1 * factor4.2 * factor4.3 * (sum(argument4.1) - sum(argument4.2))


    factor5.1 <- 1 / (8 * nu ^ 2)
    factor5.2 <- (2 * (besselI(nu, 1) / besselI(nu, 0)) ^ 2 - besselI(nu, 2) / besselI(nu, 0))
    factor5.3 <- 1 / (2 * pi * besselI(nu, 0) * n * (n - 1))

    argument5.1 <- exp(nu * cos(grid)) * (
      nu ^ 4 * sin(grid) ^ 4 - 6 * nu ^ 3 * sin(grid) ^ 2 * cos(grid) +
        3 * nu ^ 2 * (cos(grid) ^ 2 - sin(grid) ^
                        2) - nu ^ 2 * sin(grid) ^ 2  + nu * cos(grid)
    )
    argument5.2 <- exp(nu * cos(rep(0, n))) * (
      nu ^ 4 * sin(rep(0, n)) ^ 4 - 6 * nu ^ 3 * sin(rep(0, n)) ^ 2 * cos(rep(0, n)) +
        3 * nu ^ 2 * (cos(rep(0, n)) ^
                        2 - sin(rep(0, n)) ^ 2) - nu ^ 2 * sin(rep(0, n)) ^ 2  + nu * cos(rep(0, n))
    )

    part5 <- factor5.1 * factor5.2 * factor5.3 * (sum(argument5.1) - sum(argument5.2))

    output <- part2 - part3 + part4 + part5
    return(output)
  }
  bw <- optimize(
    ccv,
    interval = c(lower, upper),
    maximum = FALSE,
    x = x
  )$minimum
  if (bw < lower + tol | bw > upper - tol) {
    cli::cli_alert_warning("Minimum/maximum occurred at one end of the range")
  }
  return(bw)
}
