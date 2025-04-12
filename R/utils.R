local.factor <- function(x,
                         bw0,
                         alpha = 0.5,
                         type = c("am", "gm", "rv", "n")) {
  type <- match.arg(type)
  lambdas <- kernel.density.estimate(z = x, x = x, bw = bw0)
  g <- switch(
    type,
    "am" = mean(lambdas),
    "gm" = {
      if (any(lambdas <= 0))
        stop("For geometric mean, all density estimates must be positive.")
      exp(mean(log(lambdas)))
    },
    "rv" = max(lambdas) - min(lambdas),
    "n" = 1
  )
  return((lambdas / g) ^ (-alpha))
}

kernel.density.estimate <- function(z, x, bw) {
  n <- length(x)
  y <- 1 / (2 * n * pi * besselI(bw, 0)) * sum(exp (bw * cos(z - x)))
  return(y)
}





