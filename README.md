# pacpack

pcapack is a high-performance C library with R bindings that can
be used to quickly compute principal components, including truncated
methods.

Note that the package is currently under development, and is not
particularly stable.



## Benchmarks

Covariance:
```
     test replications elapsed relative
2 pcapack           10   0.186    1.000
1       R           10   3.456   18.581
```

PCA:
```
     test replications elapsed relative
1 pcapack           10   2.538    1.000
2       R           10   3.542    1.396
```

SVD:
```
     test replications elapsed relative
1 pcapack           10   1.550    1.000
2       R           10   1.896    1.223
```

Centering and Scaling:
```
### center=TRUE, scale=FALSE
     test replications elapsed relative
2 pcapack           10   0.244    1.000
1       R           10   2.356    9.656

### center=FALSE, scale=TRUE
     test replications elapsed relative
2 pcapack           10   0.244    1.000
1       R           10   2.356    9.656

### center=TRUE, scale=TRUE
     test replications elapsed relative
2 pcapack           10   0.497    1.000
1       R           10   7.786   15.666
```

You can find the source for these benchmarks in the `inst/benchmarks/` tree.
All tests performed using:

* R 3.2.0
* OpenBLAS
* gcc 4.9.1
* 4 cores of a Core i5-2500K CPU @ 3.30GHz


## Requirements and Installation

To install, you will need: 

* cmake >= 2.8.1
* A C99 compatible compiler with OpenMP >= 4 support.
* LAPACK and BLAS libraries (will use R's if installing the R package)
* R >= 2.14.0 and the RNACI package (if installing the R package)

Both the R package and the standalone library require cmake.

To install the R package, simply execute:

```r
library(devtools)
install_github("wrahtematics/RNACI") ### dependency
install_github("wrathematics/pcapack")
```

To build just the shared library, in your terminal, execute `make` in
`pcapack/src/pcapack/`.  A static and dynamic library will be placed 
in the `pcapack/src/pcapack/build` tree.

