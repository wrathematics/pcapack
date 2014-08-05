LA_svd <- function(x, nu=min(m, n), nv=min(m, n))
{
  m <- nrow(x)
  n <- ncol(x)
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  .Call("R_pcapack_svd", as.integer(nu), as.integer(nv), 
        as.integer(m), as.integer(n), x, 
        PACKAGE="pcapack")
  
}

