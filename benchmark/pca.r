library(rbenchmark)
library(pcapack)


set.seed(1234)
m <- 1200
n <- 350
x <- matrix(rnorm(m*n), m, n)

reps <- 2
retx <- FALSE

mine <- function(method) {mdl1 <<- pca(x, retx=retx, method=method)}
R <- function(method)
{
  if (method == "svd")
    mdl2 <<- prcomp(x, retx=retx)
  else if (method == "eigcov")
    mdl2 <<- princomp(x, retx=retx)
  
  invisible()
}

cols <- c("test", "replications", "elapsed", "relative")
benchmark(mine("svd"), R("svd"), columns=cols, replications=reps)
#all.equal(mdl1, mdl2)

benchmark(mine("eigcov"), R("eigcov"), columns=cols, replications=reps)
#all.equal(mdl1, mdl2)
