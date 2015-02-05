library(rbenchmark)
library(pcapack)


set.seed(1234)
m <- 4000
n <- 250
x <- matrix(rnorm(m*n), m, n)

reps <- 10

mine <- function() {mdl1 <<- LA_svd(x)}
R <- function() {mdl2 <<- La.svd(x)}

cols <- c("test", "replications", "elapsed", "relative")
benchmark(mine(), R(), columns=cols, replications=reps)

all.equal(mdl1, mdl2)
