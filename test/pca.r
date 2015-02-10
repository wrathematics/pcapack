library(pcapack)

set.seed(1234)

m <- 10
n <- 50
m <- 50
n <- 10

x <- matrix(rnorm(m*n), m, n)



mdl1 <- prcomp(x)
mdl2 <- pca(x, method="svd")
print(all.equal(mdl1$sdev, mdl2$sdev))
print(all.equal(mdl1$rotation, mdl2$rotation))


mdl1 <- prcomp(x)
mdl2 <- pca(x, method="eigcov")
print(all.equal(mdl1$sdev, mdl2$sdev))
print(all.equal(unclass(mdl1$loadings), mdl2$rotation))
