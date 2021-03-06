library(pcapack)
library(rbenchmark)

m <- 1e4
n <- 250
x <- matrix(rnorm(m*n), m, n)

reps <- 10

benchmark(R=A <- cov(x), pcapack=B <- cov2(x), replications=reps, columns=c("test", "replications", "elapsed", "relative"))
all.equal(A, B)
