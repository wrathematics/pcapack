library(pcapack)

set.seed(1234)
m <- 2000
n <- 500
x <- matrix(rnorm(m*n), m, n)

system.time(mdl1 <- LA_svd(x))[3]
system.time(mdl2 <- La.svd(x))[3]

all.equal(mdl1, mdl2)
