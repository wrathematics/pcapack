#include <stdlib.h>
#include "benchmark.h"
#include "../src/pcapack.h"

#define N     6000
#define REPS  5

int main()
{
  double *x = malloc(N*N * sizeof(double));
  
  
  BENCHMARK(pcapack_xpose(N, N, x), REPS);
  
  free(x);
  return 0;
}
