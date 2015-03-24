library(pcapack)
library(rbenchmark)

set.seed(1234)

m <- 1e4
n <- 250
mean <- 0
sd <- 1
x <- matrix(rnorm(m*n, sd=sd, mean=mean), m, n)

reps <- 5


print("center=TRUE, scale=FALSE")
benchmark(R=A <- scale(x, TRUE, FALSE), Me=B <- scale2(x, TRUE, FALSE), replications=reps, columns=c("test", "elapsed", "relative"))
all.equal(A, B, check.attributes=FALSE)

print("center=FALSE, scale=TRUE")
benchmark(R=A <- scale(x, FALSE, TRUE), Me=B <- scale2(x, FALSE, TRUE), replications=reps, columns=c("test", "elapsed", "relative"))
all.equal(A, B, check.attributes=FALSE)

print("center=TRUE, scale=TRUE")
benchmark(R=A <- scale(x, TRUE, TRUE), Me=B <- scale2(x, TRUE, TRUE), replications=reps, columns=c("test", "elapsed", "relative"))
all.equal(A, B, check.attributes=FALSE)
