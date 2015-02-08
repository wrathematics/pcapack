library(pcapack)

m <- 100
n <- 25
x <- matrix(rnorm(m*n), m, n)

all.equal(cov(x), pcapack_cov(x))

