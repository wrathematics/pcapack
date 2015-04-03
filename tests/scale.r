library(pcapack)

set.seed(1234)

m <- 100
n <- 25
mean <- 0
sd <- 1
x <- matrix(rnorm(m*n, sd=sd, mean=mean), m, n)


# center=TRUE, scale=FALSE
A <- scale(x, TRUE, FALSE)
B <- scale2(x, TRUE, FALSE)
stopifnot(all.equal(A, B, check.attributes=FALSE))

# center=FALSE, scale=TRUE
A <- scale(x, FALSE, TRUE)
B <- scale2(x, FALSE, TRUE)
stopifnot(all.equal(A, B, check.attributes=FALSE))

# center=TRUE, scale=TRUE
A <- scale(x, TRUE, TRUE)
B <- scale2(x, TRUE, TRUE)
stopifnot(all.equal(A, B, check.attributes=FALSE))

