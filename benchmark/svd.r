library(rbenchmark)
library(pcapack)


set.seed(1234)
m <- 5000
n <- 250
x <- matrix(rnorm(m*n), m, n)

reps <- 5

mine <- function() {mdl1 <<- LA_svd(x)}
R <- function() {mdl2 <<- La.svd(x)}

cols <- c("test", "replications", "elapsed", "relative")
benchmark(Me = mine(), R = R(), columns=cols, replications=reps)

all.equal(mdl1, mdl2)
