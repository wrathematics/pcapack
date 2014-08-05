LA_svd <- function(x, nu=min(m, n), nv=min(m, n))
{
  m <- nrow(x)
  n <- ncol(x)
  k <- min(m, n)
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  ret <- .Call("R_pcapack_svd", as.integer(nu), as.integer(nv), 
        as.integer(m), as.integer(n), x, 
        PACKAGE="pcapack")
  
  if (nu)
  {
    if (nu != k && nu != m)
      ret$u <- ret$u[, 1L:nu]
  }
  
  if (nv)
  {
    if (nv != k && nv != n)
      ret$vt <- ret$vt[1L:nv, ]
  }
  
  return( ret )
}

