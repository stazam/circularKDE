---
title: "README"
author: "Stanislav Zámečník"
date: "15 04 2025"
output:
  html_document: default
  pdf_document: default
---
# circularKDE

## Methods for kernel density estimation for circular data

The original version of this software was written in R by Stanislav Zámečník with the help of Stanislav Katina in year 2025. The package is based on theoretical background of work of prof. Ivanka Hórová, Kamila Hasilová, Stanislav Zámečník and afterwards implemented by mentioned authors. 

Main features of the package include:
- new methods implementation (see DESCRIPTION).
- adaptive density estimation function


## Installation

You can install the release version of package from [CRAN](https://CRAN.R-project.org): 

``` r
install.packages("circularKDE")
``` 

Or the development version from GitHub repository:

``` r
devtools::install_github("stazam/circularKDE")
``` 

## Usage

In R session do:

```r
library(circular)
x_circ <- rvonmises(100, mu = circular(0), kappa = 1)
```

To use estimation by adaptive density function simply run 

```r
library(circularKDE)

bw <- bw.lscvg(x_circ)
dens <- density.adaptive.circular(x, bw0 = bw0)
plot(seq(0, 2 * pi, length.out = 500), dens, type = "l",
      main = "Adaptive Circular Density")
```


For more information and examples see:

```r
?circularKDE
```
