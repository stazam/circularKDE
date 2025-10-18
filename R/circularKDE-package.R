#' @title circularKDE: Recent Methods for Kernel Density Estimation of Circular Data
#'
#' @description
#' The circularKDE package provides state-of-the-art nonparametric kernel density 
#' estimation methods specifically designed for circular (directional) data. The package 
#' implements recent advances in bandwidth selection algorithms and adaptive density 
#' estimation techniques for analyzing angular and directional datasets.
#'
#' @details
#' The package includes several bandwidth selection methods:
#' \itemize{
#'   \item \code{\link{bwScv}} - Smoothed Cross-Validation bandwidth selection
#'   \item \code{\link{bwCcv}} - Complete Cross-Validation method  
#'   \item \code{\link{bwLscvg}} - Generalized Least Squares Cross-Validation
#'   \item \code{\link{bwFo}} - Fourier-based plug-in bandwidth selection
#'   \item \code{\link{bwJf}} - Jones-Foster additive plug-in method
#'   \item \code{\link{bwTs}} - Terrell-Scott multiplicative plug-in approach
#' }
#'
#' The main density estimation function is:
#' \itemize{
#'   \item \code{\link{adaptiveDensityCircular}} - Variable bandwidth kernel density estimation
#' }
#'
#' All methods are optimized for computational efficiency while maintaining numerical 
#' accuracy for both small and large datasets. The implementation is based on recent 
#' theoretical advances by Hórová, Hasilová, Zámečník, and Katina (2024-2025).
#'
#' @section Getting Started:
#' To get started with circularKDE:
#' \preformatted{
#' library(circular)
#' library(circularKDE)
#' 
#' # Generate some circular data
#' x_circ <- rvonmises(100, mu = circular(0), kappa = 1)
#' 
#' # Select bandwidth using LSCV-g method
#' bw <- bwLscvg(x_circ)
#' 
#' # Estimate adaptive density
#' dens <- adaptiveDensityCircular(x_circ, bw0 = bw)
#' 
#' # Plot the result
#' plot(seq(0, 2 * pi, length.out = 500), dens, type = "l")
#' }
#'
#' @author 
#' Stanislav Zámečník \email{zamecnik@@math.muni.cz}
#' 
#' Ivanka Hórová \email{horova@@math.muni.cz}
#' 
#' Kamila Hasilová \email{kamila.hasilova@@unob.cz}
#' 
#' Stanislav Katina \email{katina@@math.muni.cz}
#'
#' @references
#' Zámečník, S., Horová, I., Katina, S., & Hasilová, K. (2024). 
#' An adaptive method for bandwidth selection in circular kernel density estimation. 
#' \emph{Computational Statistics}. \doi{10.1007/s00180-023-01401-0}
#' 
#' Hasilová, K., Horová, I., Valis, D., & Zámečník, S. (2024). 
#' A comprehensive exploration of complete cross-validation for circular data. 
#' \emph{Statistics in Transition New Series}, 25(3):1--12. \doi{10.59170/stattrans-2024-024}
#'
#' @keywords package
#' @docType package
#' @name circularKDE-package
#' @aliases circularKDE circularKDE-package
#' @importFrom circular circular
#' @importFrom cli cli_abort cli_warn
NULL
