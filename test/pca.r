library(pcapack)

set.seed(1234)

m <- 10
n <- 5
k <- min(m, n)

x <- matrix(rnorm(m*n), m, n)


test <- function()
{
  mdl1 <- prcomp(x)
  mdl2 <- pca(x)
  
#  all.equal(mdl1, mdl2)
  print(all.equal(mdl1$sdev, mdl2$sdev))
  print(all.equal(mdl1$rotation, mdl2$rotation))
  
  invisible()
}



test()

