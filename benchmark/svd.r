library(rbenchmark)
library(pcapack)


set.seed(1234)
m <- 5000
n <- 350
x <- matrix(rnorm(m*n), m, n)

reps <- 10

mine <- function() {mdl1 <<- svd2(x)}
R <- function() {mdl2 <<- La.svd(x)}

cols <- c("test", "replications", "elapsed", "relative")
benchmark(pcapack=mine(), R=R(), columns=cols, replications=reps)

all.equal(mdl1, mdl2)
