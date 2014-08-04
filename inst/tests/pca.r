library(ffpca)

n <- 3e4
p <- 1e2

x <- matrix(rnorm(n*p), n, p)


system.time(prcomp(x))[3]
system.time(pca(x))[3]
