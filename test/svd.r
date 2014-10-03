library(pcapack)

set.seed(1234)

test <- function(x, nu, nv)
{
  m1 <- LA_svd(x, nu, nv)
  m2 <- La.svd(x, nu, nv)
  
  check <- all.equal(m1, m2)
  
  cat(paste("nu=", nu, " nv=", nv, ":\t", check, "\n", sep=""))
  invisible()
}

suite <- function(m, n)
{
  k <- min(m, n)
  x <- matrix(rnorm(m*n), m, n)
  
  test(x, 0, 0)
  test(x, k, k)
  test(x, k, 0)
  test(x, 0, k)
  test(x, m, 0)
  test(x, 0, n)
  test(x, m, n)
  test(x, 5, 3)
  
  invisible()
}

#########################


suite(10, 5)
#suite(5, 10)
