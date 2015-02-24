#include <stdlib.h>
#include "benchmark.h"
#include "../src/pcapack.h"

#define M     20000
#define N     250
#define REPS  15

int main()
{
  double *x = malloc(M*N * sizeof(double));
  double *cov = malloc(N*N * sizeof(double));
  
  BENCHMARK(pcapack_cov(COR_PEARSON, M, N, x, cov), REPS);
  
  free(x);
  free(cov);
  return 0;
}
