#ifndef PCAPACK_RAND_H
#define PCAPACK_RAND_H


#define I32BIT 2.328306437080797e-10
#define MT_SIZE 624
#define ZT_SIZE 256

void mt_check();
void mt_init_gen(unsigned int seed);
void mt_gen();
unsigned int mt_extract();

int sample1(const int min, const int max);
void sample1_(int* min, int* max, int* ret);
double runif1(const double min, const double max);
void rnormn_(int *n, double *mn, double *sd, double *ret);


#endif
