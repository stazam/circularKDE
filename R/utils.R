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
