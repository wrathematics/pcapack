library(rbenchmark)
library(pcapack)


set.seed(1234)
m <- 1200
n <- 350
x <- matrix(rnorm(m*n), m, n)

reps <- 2

mine <- function() {mdl1 <<- pca(x)}
R <- function() {mdl2 <<- prcomp(x)}

cols <- c("test", "replications", "elapsed", "relative")
benchmark(mine(), R(), columns=cols, replications=reps)


#all.equal(mdl1, mdl2)
