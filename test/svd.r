library(pcapack)

set.seed(1234)

m <- 10
n <- 5
k <- min(m, n)

x <- matrix(rnorm(m*n), m, n)


test <- function(nu, nv)
{
  LA_svd(x, nu, nv)
  La.svd(x, nu, nv)
  
  all.equal(LA_svd(x), La.svd(x))
}



test(0, 0)
test(k, k)
test(k, 0)
test(0, k)
test(m, 0)
test(0, n)
test(m, n)
