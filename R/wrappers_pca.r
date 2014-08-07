pca <- function(x, retx=TRUE, center=TRUE, scale.=FALSE)
{
  m <- nrow(x)
  n <- ncol(x)
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  .Call("R_pca", as.integer(m), as.integer(n), x, 
        as.integer(center), as.integer(scale.), as.integer(retx),
        PACKAGE="pcapack")
  
}

