pcapack_scale <- function(x, center=TRUE, scale=TRUE)
{
  if (!is.matrix(x))
    dim(x) <- c(length(x), 1L)
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  .Call(R_scale, as.logical(center), as.logical(scale), x)
}
