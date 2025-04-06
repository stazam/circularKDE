#' Title
#'
#' @param x something
#' @param bw0 something
#' @param alpha something
#' @param type something
#' @param z something
#' @param from something
#' @param to something
#' @param n something
#'
#' @return something
#' @export
#'
#' @examples something
density.adaptive.circular <- function(x,
                                      bw0,
                                      alpha = 0.5,
                                      type = "n",
                                      z = NULL,
                                      from = circular(0),
                                      to = circular(2 * pi),
                                      n = 500) {
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
  if (!is.numeric(from)) {
    cli_abort("{.var from} must be a numeric argument")
  }
  if (!is.finite(from)) {
    cli_abort("{.var from} must be a non-fininte argument")
  }
  if (!is.numeric(to)) {
    cli_abort("{.var from} must be a numeric argument")
  }
  if (!is.finite(to)) {
    cli_abort("{.var from} must be a non-fininte argument")
  }
  n <- round(n)
  if (!is.numeric(n)) {
    cli_abort("{.var from} must be a numeric argument")
  }
  if (n <= 0) {
    cli_abort("argument 'n' must be integer and positive")
  }
  if (is.null(z)) {
    z <- circular(seq(
      from = from,
      to = to,
      length = n
    ))
  }
  density.est <- function(z, x, bw0, alpha, type) {
    n <- length(x)
    factor <- 1 / (2 * n * pi)

    lambda <- local.factor(x, bw0, alpha, type)
    main.part <- sum(1 / besselI(lambda * bw0, 0) * exp (lambda * bw0 * cos(z - x)))
    return(factor * main.part)
  }
  y <- sapply(z,
              density.est,
              x = x,
              bw0 = bw0,
              alpha = alpha,
              type = type)
  return(y)
}

#' @export
density.hat <- function(z, x, bw) {
  n <- length(x)
  y <- 1 / (2 * n * pi * besselI(bw, 0)) * sum(exp (bw * cos(z - x)))
  return(y)
}

#' @export
local.factor <- function(x,
                         bw0,
                         alpha = 0.5,
                         type = c("am", "gm", "rv", "n")) {
  lambdas <- sapply(x, density.hat, x = x, bw = bw0)
  g <- switch(
    type,
    "am" = mean(lambdas),
    "gm" = exp(mean(log(lambdas))),
    "rv" = max(lambdas) - min(lambdas),
    "n" = 1
  )
  return((lambdas / g) ^ (-alpha))

}
