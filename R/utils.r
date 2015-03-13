is.int <- function(x)
{
  if (is.numeric(x))
  {
    if (x-as.integer(x) == 0)
      return( TRUE )
    else
      return( FALSE )
  }
  else
    return( FALSE )
}



check_mvdf <- function(x)
{
  if (is.data.frame(x))
  {
    expr <- substitute(x <- as.matrix(x))
    eval(expr, parent.frame())
  }
  else if (is.numeric(x) && !is.matrix(x))
    stop("input must be a matrix, vector, or dataframe")
}

