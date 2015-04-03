#' SVD
#' 
#' Singular value decomposition.
#' 
#' @param x
#' The input matrix or dataframe.
#' @param nu
#' Number of left singular vectors to return.
#' @param nv
#' Number of right singular vectors to return.
#' 
#' @rdname svd
#' @export
svd2 <- function(x, nu=min(n, p), nv=min(n, p))
{
  assert.type(x, "numeric")
  assert.natnum(nu)
  assert.natnum(nv)

  if (!is.matrix(x))
    x <- as.matrix(x)
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  n <- nrow(x)
  p <- ncol(x)
  
  ret <- .Call(R_pcapack_svd, as.integer(nu), as.integer(nv), x)
  
  k <- min(dim(x))
  
  if (nu)
  {
    if (nu != k && nu != n)
      ret$u <- ret$u[, 1L:nu]
  }
  
  if (nv)
  {
    if (nv != k && nv != p)
      ret$vt <- ret$vt[1L:nv, ]
  }
  
  return( ret )
}

