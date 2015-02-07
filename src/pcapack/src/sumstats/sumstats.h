#ifndef __PCAPACK_SUMSTATS_H__
#define __PCAPACK_SUMSTATS_H__

#include <stdbool.h>

// crossprod.h
#define UPPER 1
#define LOWER 2

int pcapack_symmetrize(const int triang, const int m, const int n, double *x);
int pcapack_crossprod(int m, int n, double *x, double alpha, double *c);
int pcapack_tcrossprod(int m, int n, double *x, double alpha, double *c);
int pcapack_inverse(int n, double *x);



// covariance.c
#define COR_PEARSON   1
#define COR_SPEARMAN  2
#define COR_KENDALL   3

int pcapack_cov(const int method, int m, int n, double *x, double *cov);
int pcapack_cor(const int method, int m, int n, double *x, double *cor);



// sweeps.c
#define PLUS 1
#define MINUS 2
#define TIMES 3
#define DIVIDE 4

#define ROWS 1
#define COLS 2

int pcapack_sweep(const int m, const int n, double *x, double *vec, int lvec, int margin, int fun);
int pcapack_scale(bool centerx, bool scalex, const int m, const int n, double *x);


#endif
