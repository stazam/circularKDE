localFactor <- function(x,
                         bw0,
                         alpha = 0.5,
                         type = c("am", "gm", "rv", "n")) {
  type <- match.arg(type)
  lambdas <- vapply(x, function(z_i) kernelDensityEst(z_i, x = x, bw = bw0), numeric(1))
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
  result <- (g / lambdas) ^ (alpha)
  return(result)
}

kernelDensityEst <- function(z, x, bw) {
  n <- length(x)
  factor <- 1 / (2 * n * pi * besselI(bw, 0))
  result <- factor * rowSums(exp(bw * cos(outer(z, x, "-"))))
  return(result)
}



