# pacpack

pcapack is a combination shared library and R package that can
be used to quickly compute principal components, including using
approximate truncated methods.



## Requirements and Installation

To install, you will need: 

* cmake >= 2.8.1
* A Fortran 2003 compatible compiler (for iso_c_binding) with OpenMP support.
* A C99 compatible compiler
* LAPACK and BLAS libraries (will use R's if installing the R package)
* R >= 2.14.0 and the RNACI package (if installing the R package)

If your make -j defaults to something greater than 1, the compile
may break; try again 1 or 2 times, or use make -j 1.

Both the R package and the standalone library require cmake, because if you
so much as think the word "autotools" around me, I'll punch you in the 
stomach.

To install the R package, simply execute:

```
R CMD INSTALL pcapack_0.1.0.tar.gz
```

To build just the shared library, in your terminal, navigate to:

```
cd pcapack/src/pcapack/
make
```

A static and dynamic library will be placed in the 
`pcapack/src/pcapack/build` tree.


