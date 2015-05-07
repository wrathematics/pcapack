library(pcapack)

set.seed(1234)

m <- 10
n <- 50
m <- 50
n <- 10

x <- matrix(rnorm(m*n), m, n)



mdl1 <- prcomp(x)
mdl1$center <- mdl1$scale <- NULL
mdl2 <- pca(x, method="svd")
stopifnot(all.equal(mdl1, mdl2))



#princomp2prcomp <- function(mdl)
#{
  #names(mdl1$sdev) <- NULL
  #rownames(mdl1$loadings) <- NULL
  #colnames(mdl1$loadings) <- 
  
#}


#mdl1 <- princomp(x)
#mdl1 <- princomp2prcomp(mdl1)
#mdl2 <- pca(x, method="eigcov")
#stopifnot(all.equal(mdl1, mdl2))

#print(mdl2)
