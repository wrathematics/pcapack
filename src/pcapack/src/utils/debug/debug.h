#ifndef __PCAPACK_DEBUG_H__
#define __PCAPACK_DEBUG_H__


#include <stdbool.h>

// compare_mats.c
bool isequal_mats(int m, int n, double *x, double *y, double tol);

// matprinter.c
void matprinter(int m, int n, double *x);


#endif
