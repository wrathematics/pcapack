#' PCA
#' 
#' Principal components analysis.
#' 
#' @param x
#' The input matrix or dataframe.
#' @param retx
#' Return the rotated variables?
#' @param center
#' Center the matrix first?
#' @param scale
#' Scale the matrix first?
#' @param method
#' "svd" for svd of data matrix, or "eigcov" for eigenvalue decomposition
#' of the covariance matrix.
#' 
#' @rdname pca
#' @export
pca <- function(x, retx=TRUE, center=TRUE, scale=FALSE, method="svd")
{
  assert.type(method, "character")
  method <- match.arg(tolower(method), c("svd", "eigcov"))
  
  assert.type(retx, "logical")
  assert.type(center, "logical")
  assert.type(scale, "logical")
  
  assert.type(x, "numeric")

  if (!is.double(x))
    storage.mode(x) <- "double"
  
  if (method == "svd")
    ret <- .Call(R_pcapack_prcomp_svd, x, as.integer(center), as.integer(scale), as.integer(retx))
  else if (method == "eigcov")
    ret <- .Call(R_pcapack_prcomp_eigcov, x, as.integer(retx))
  
  class(ret) <- "prcomp"
  
  return(ret)
}

