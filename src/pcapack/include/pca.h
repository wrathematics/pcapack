#ifndef PCAPACK_PCA_H
#define PCAPACK_PCA_H


void dpca_(int *m, int *n, int *k, double *x, double *sdev, double *trot, char *ret, char *cntr, char *scl, int *info);
#define DPCA dpca_


#endif
