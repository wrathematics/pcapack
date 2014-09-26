cma <- function(x, k=1)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  .Call("R_cma", x, as.integer(k), PACKAGE="pcapack")
}

