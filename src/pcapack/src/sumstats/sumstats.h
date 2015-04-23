#ifndef __PCAPACK_SUMSTATS_H__
#define __PCAPACK_SUMSTATS_H__

#include <stdbool.h>


// covariance.c
#define COR_PEARSON   1
#define COR_SPEARMAN  2
#define COR_KENDALL   3

int pcapack_cov(const int method, int m, int n, double *x, double *cov);
int pcapack_cor(const int method, int m, int n, double *x, double *cor);



// crossprod.c
#define UPPER 1
#define LOWER 2

bool pcapack_is_symmetric(const int m, const int n, const double *x);
int pcapack_symmetrize(const int triang, const int m, const int n, double *x);
int pcapack_crossprod(int m, int n, const double *x, double alpha, double *c);
int pcapack_tcrossprod(int m, int n, const double *x, double alpha, double *c);



// means.c
int pcapack_rowsums(const int m, const int n, double *x, double *rowsums);
int pcapack_colsums(const int m, const int n, double *x, double *colsums);
double pcapack_mean(const int n, double *x);
int pcapack_rowmeans(const int m, const int n, double *x, double *rowsums);
int pcapack_colmeans(const int m, const int n, double *x, double *colsums);



// sweeps.c
#define PLUS 1
#define MINUS 2
#define TIMES 3
#define DIVIDE 4

#define ROWS 1
#define COLS 2

int pcapack_sweep(const int m, const int n, double *x, double *vec, int lvec, int margin, int fun);
int pcapack_scale(const bool centerx, const bool scalex, const int m, const int n, double *x);



// variances.c
double pcapack_variance(const int n, double *x);
double pcapack_sdev(const int n, double *x);
int pcapack_rowvars(const int m, const int n, double *x, double *rowvars);
int pcapack_colvars(const int m, const int n, double *x, double *colvars);
int pcapack_colsdev(const int m, const int n, double *x, double *colsdev);


#endif
