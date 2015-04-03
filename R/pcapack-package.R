#' PCAPACK
#' 
#' PCAPACK is a set of utilities for performing principal components
#' analysis (including by means of truncated svd) very efficiently.
#' The package also contains efficient, related helper utilities
#' such as a replacement for R's covariance calculator.
#' 
#' \tabular{ll}{ 
#'    Package: \tab PCAPACK \cr 
#'    Type: \tab Package \cr 
#'    License: \tab MPL \cr 
#'    Needs Compilation: \tab yes \cr
#' } 
#' 
#' 
#' @useDynLib pcapack, R_cov, R_scale, R_pcapack_prcomp_svd, 
#' R_pcapack_prcomp_eigcov, R_pcapack_svd
#' 
#' @name PCAPACK-package
#' @docType package
#' @author Drew Schmidt
#' @keywords Package
NULL

