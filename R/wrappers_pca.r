pca <- function(x, retx=TRUE, center=TRUE, scale=FALSE, method="svd")
{
  method <- match.arg(tolower(method), c("svd", "eigcov"))
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  if (method == "svd")
    ret <- .Call("R_pca_svd", x, as.integer(center), as.integer(scale), as.integer(retx), PACKAGE="pcapack")
  else if (method == "eigcov")
    ret <- .Call("R_pca_eigcov", x, as.integer(retx), PACKAGE="pcapack")
  
  return(ret)
}

