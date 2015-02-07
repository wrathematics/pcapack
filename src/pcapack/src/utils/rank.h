#ifndef __PCAPACK_RANK_H__
#define __PCAPACK_RANK_H__


#define RANK_MIN 1
#define RANK_MAX 2
#define RANK_AVG 3
#define RANK_DEF 4

void genrank(const double *data, double *rank, const int len, int method);


#endif
