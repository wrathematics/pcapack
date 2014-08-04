cma <- function(x, k=1)
{
  m <- nrow(x)
  n <- ncol(x)
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  .Call("R_cma",
        as.integer(m), as.integer(n), x, as.integer(k),
        PACKAGE = "ffpca")
}

