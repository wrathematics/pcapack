pca <- function(x, center=TRUE, scale=FALSE, retrot=TRUE)
{
  m <- nrow(x)
  n <- ncol(x)
  k <- min(m, n)
  
  if (center)
    center <- 'Y'
  else
    center <- 'N'
  
  if (scale)
    scale <- 'Y'
  else
    scale <- 'N'
  
  if (retrot)
    retrot <- 'Y'
  else
    retrot <- 'N'
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  .Call("R_pca",
        as.integer(m), as.integer(n), as.integer(k), x, as.character(center), as.character(scale), as.character(retrot),
        PACKAGE = "ffpca")
}

