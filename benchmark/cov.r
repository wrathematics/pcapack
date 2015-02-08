library(pcapack)
library(rbenchmark)

m <- 1e4
n <- 250
x <- matrix(rnorm(m*n), m, n)

reps <- 5

benchmark(cov(x), pcapack_cov(x), replications=reps, columns=c("test", "elapsed", "relative"))

