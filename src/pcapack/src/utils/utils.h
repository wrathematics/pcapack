#ifndef __PCAPACK_UTILS_H__
#define __PCAPACK_UTILS_H__

#include "rank.h"


// distances.c
double pcapack_distance(const int method, const int n, const double *x, const double *y, const int p);
double pcapack_mahalanobis(int n, const double *x, const double *y, const double *covinv);


// dnrm2.c
double dnrm2(int n, double *x, int incx);


// inverse.c
int pcapack_inverse(int n, double *x);


// Rspecial.c
bool pcapack_anyna(const int n, const double *x);


// xpose.c
void pcapack_xpose(int m, int n, double *x);


#endif


