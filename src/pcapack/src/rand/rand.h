#ifndef PCAPACK_RAND_H
#define PCAPACK_RAND_H


#include "rng_interface.h"

#define I32BIT 2.328306437080797e-10
#define MT_SIZE 624
#define ZT_SIZE 256


// rand_normal.c
double stdnorm(rng_state_t *rs);
double rnorm(rng_state_t *rs, const double mn, const double sd);


// rand_uniform.c
double runif(rng_state_t *rs, const double min, const double max);


// samplers.c
int sample(rng_state_t *rs, const int min, const int max);


#endif

