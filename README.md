# pacpack

pcapack is a high-performance C library with R bindings that can
be used to quickly compute principal components, including truncated
methods.



## Requirements and Installation

To install, you will need: 

* cmake >= 2.8.1
* A C99 compatible compiler with OpenMP support.
* LAPACK and BLAS libraries (will use R's if installing the R package)
* R >= 2.14.0 and the RNACI package (if installing the R package)

Both the R package and the standalone library require cmake, because if you
so much as think the word "autotools" around me, I'll punch you in the 
stomach.

To install the R package, simply execute:

```
library(devtools)
install_github("wrahtematics/RNACI") ### dependency
install_github("wrathematics/pcapack")
```

To build just the shared library, in your terminal, execute `make` in
`pcapack/src/pcapack/`.  A static and dynamic library will be placed 
in the `pcapack/src/pcapack/build` tree.

