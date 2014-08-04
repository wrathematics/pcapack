# Copyright 2013, Schmidt

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

