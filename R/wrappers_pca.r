pca <- function(x, center=TRUE, scale=FALSE, retrot=TRUE)
{
  m <- nrow(x)
  n <- ncol(x)
  k <- min(m, n)
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  .Call("R_pca", as.integer(m), as.integer(n), as.integer(k), x, 
        as.integer(center), as.integer(scale), as.integer(retrot),
        PACKAGE="pcapack")
  
}

