#ifndef __PCAPACK_UTILS_H__
#define __PCAPACK_UTILS_H__

#include "rank.h"

// dnrm2.c
double dnrm2(int n, double *x, int incx);

// matprinter.c
void matprinter(int m, int n, double *x);

// Rspecial.c
bool pcapack_anyna(const int n, const double *x);

// xpose.c
void pcapack_xpose(int m, int n, double *x);


#endif


