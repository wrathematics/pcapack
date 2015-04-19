#' CMA
#' 
#' Coordinate Mapping Algorithm.
#' 
#' @param x
#' The input matrix or dataframe.
#' @param k
#' TODO (number of subspace iterations I think, don't remember)
#' 
#' @rdname cma
#' @export
cma <- function(x, k=1)
{
  assert.natnum(k)
  assert.type(x, "numeric")
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  ret <- .Call(R_pcapack_cma, x, as.integer(k))
  
  return(ret)
}

