#' Covariance
#' 
#' An efficient method of quickly computing the covariance of
#' a matrix.  Like R's own \code{cov()}, this will cast dataframes 
#' as matrices first, which is potentially a very memory expensive
#' operation.
#' 
#' @param x
#' The input matrix, vector, or dataframe.
#' 
#' @rdname cov
#' @export
cov2 <- function(x)
{
  check_mvdf(x)
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  .Call(R_cov, x)
}
