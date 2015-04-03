#' Scale
#' 
#' An efficient method of centering and/or scaling a matrix or vector.
#' Like R's own \code{scale()}, this will cast dataframes as matrices
#' first, which is potentially a very memory expensive operation.
#' 
#' @param x
#' The input matrix, vector, or dataframe.
#' 
#' @rdname scale
#' @export
scale2 <- function(x, center=TRUE, scale=TRUE)
{
  assert.type(center, "logical")
  assert.type(scale, "logical")
  check_mvdf(x)
  
  assert.type(x, "numeric")
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  .Call(R_scale, as.logical(center), as.logical(scale), x)
}

