# circularKDE 
<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/circularKDE)](https://CRAN.R-project.org/package=circularKDE)
[![R-CMD-check](https://github.com/stazam/circularKDE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/stazam/circularKDE/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/stazam/circularKDE/branch/main/graph/badge.svg)](https://app.codecov.io/gh/stazam/circularKDE?branch=main)
[![License: GPL v2](https://img.shields.io/badge/License-GPL_v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![Downloads](https://cranlogs.r-pkg.org/badges/circularKDE)](https://CRAN.R-project.org/package=circularKDE)
<!-- badges: end -->

## Overview

Package provides state-of-the-art nonparametric kernel density estimation methods specifically designed for circular (directional) data. The package implements recent advances in bandwidth selection algorithms and adaptive density estimation techniques. Namely:

- **`bwScv()`** - Smoothed Cross-Validation bandwidth selection from Zamecnik et al. (2024)
- **`bwCcv()`** - Complete Cross-Validation method derived from Hasilová et al. (2024) for circular kernel density estimation  
- **`bwLscvg()`** - Generalized Least Squares Cross-Validation from Zamecnik et al. (2025)
- **`bwFo()`** - Fourier-based plug-in bandwidth selection following Tenreiro (2022) methodology
- **`bwJf()`** - Jones-Foster additive plug-in method adapted from Tsuruta et al. (2017) for circular data
- **`bwTs()`** - Terrell-Scott multiplicative plug-in approach based on Tsuruta et al. (2017) circular modifications
- **`adaptiveDensityCircular()`** - Variable bandwidth kernel density estimation using local adaptation factors

Additional information on the choice of technique and comparison of fixed vs adaptive methods for different datasets please refer to a thesis Zamecnik (2025).



## Installation

You can install the release version of package from [CRAN](https://CRAN.R-project.org): 

``` r
install.packages("circularKDE")
``` 

Or the development version from GitHub repository [link](https://github.com/stazam/circularKDE):

``` r
devtools::install_github("stazam/circularKDE")
``` 

## Usage

### Basic Example

In R session do:

```r
library(circular)
library(circularKDE)

# sample from Von-Mises Fisher distribution
x_circ <- rvonmises(100, mu = circular(0), kappa = 1)

# computation of the bandwidth based on LSCVg method
bw <- bwLscvg(x_circ)

# evaluation and plot of the adaptive kernel density estiamtor 
dens <- adaptiveDensityCircular(x_circ, bw0 = bw)
plot(seq(0, 2 * pi, length.out = 500), dens, type = "l",
      main = "Adaptive Circular Density")
```

### Bandwidth Selection Comparison

```r
library(circular)
library(circularKDE)

set.seed(42)
x_circular <- rvonmises(100, mu = circular(0), kappa = 1)

# Compare usage of different bandwidth selection methods
(h_SCV <- bwScv(x_circular)$minimum)
(h_CCV <- bwCcv(x_circular)$minimum)
(h_LSCVg <- bwLscvg(x_circular))
(h_FO <- bwFo(x_circular))
(h_JF <- bwJf(x_circular))
(h_TS <- bwTs(x_circular))

theta <- seq(0, 2 * pi, length.out = 500)
y_SCV <- density(x = x_circular, z = theta, bw=h_SCV)$y
y_CCV <- density(x = x_circular, z = theta, bw=h_CCV)$y
y_LSCVg <- density(x = x_circular, z = theta, bw=h_LSCVg)$y
y_FO <- density(x = x_circular, z = theta, bw=h_FO)$y
y_JF <- density(x = x_circular, z = theta, bw=h_JF)$y
y_TS <- density(x = x_circular, z = theta, bw=h_TS)$y

# plot the comparison on the interval [0, 2pi]
plot(theta, y_SCV, type = "l", lwd = 2, col = "black",
     xlim = c(0, 2*pi), ylim = c(0,1),
     xlab = expression(theta), ylab = "Density",
     main = "Comparison of Circular KDE Bandwidth Selectors")

lines(theta, y_CCV,   col = "red",   lwd = 2, lty = 2)
lines(theta, y_LSCVg, col = "blue",  lwd = 2, lty = 3)
lines(theta, y_FO,    col = "green", lwd = 2, lty = 4)
lines(theta, y_JF,    col = "purple",lwd = 2, lty = 5)
lines(theta, y_TS,    col = "orange",lwd = 2, lty = 6)

legend("topright",
       legend = c("SCV", "CCV", "LSCVg", "FO", "JF", "TS"),
       col = c("black", "red", "blue", "green", "purple", "orange"),
       lwd = 2, lty = 1:6)
```

For more information and examples see:

```r
?circularKDE
```

## References

Hasilová, K., Horová, I., Valis, D., & Zámečník, S. (2024). A comprehensive exploration of complete cross-validation for circular data. *Statistics in Transition New Series*, 25(3):1--12. [doi:10.59170/stattrans-2024-024](https://doi.org/doi:10.59170/stattrans-2024-024)

Zámečník, S., Horová, I., Katina, S., & Hasilová, K. (2024). An adaptive method for bandwidth selection in circular kernel density estimation. *Computational Statistics*. 
[doi:10.1007/s00180-023-01401-0](https://doi.org/10.1007/s00180-023-01401-0)

Tenreiro, C. (2022). Kernel density estimation for circular data: a Fourier series-based plug-in approach for bandwidth selection. *Journal of Nonparametric Statistics*, 34(2):377--406. [doi:10.1080/10485252.2022.2057974](https://doi.org/10.1080/10485252.2022.2057974)

Tsuruta, Y., & Sagae, M. (2017). Higher order kernel density estimation on the circle. *Statistics & Probability Letters*, 131:46--50. [doi:10.1016/j.spl.2017.08.003](https://doi.org/10.1016/j.spl.2017.08.003)

Jones, M. C., & Foster, P. J. (1993). Generalized jackknifing and higher-order kernels. *Journal of Nonparametric Statistics*, 3:81--94. [doi:10.1080/10485259308832573](https://doi.org/10.1080/10485259308832573)

Terrell, G. R., & Scott, D. W. (1980). On improving convergence rates for nonnegative kernel density estimators. *The Annals of Statistics*, 8(5):1160--1163.




