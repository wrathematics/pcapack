pcapack_cov <- function(x)
{
  if (!is.matrix(x))
    dim(x) <- c(length(x), 1L)
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  .Call(R_cov, x)
}
