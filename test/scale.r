library(pcapack)

set.seed(1234)

m <- 100
n <- 25
mean <- 0
sd <- 1
x <- matrix(rnorm(m*n, sd=sd, mean=mean), m, n)


cat("center=TRUE, scale=FALSE\n")
A <- scale(x, TRUE, FALSE)
B <- scale2(x, TRUE, FALSE)
all.equal(A, B, check.attributes=FALSE)

cat("center=FALSE, scale=TRUE\n")
A <- scale(x, FALSE, TRUE)
B <- scale2(x, FALSE, TRUE)
all.equal(A, B, check.attributes=FALSE)

cat("center=TRUE, scale=TRUE\n")
A <- scale(x, TRUE, TRUE)
B <- scale2(x, TRUE, TRUE)
all.equal(A, B, check.attributes=FALSE)
