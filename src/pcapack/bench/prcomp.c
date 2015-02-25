#include <stdlib.h>
#include "benchmark.h"
#include "../src/pcapack.h"

#define M     5000
#define N     250
#define REPS  15

int main()
{
  const int k = M<N?M:N;
  double *x = malloc(M*N * sizeof(double));
  double *sdev = malloc(k * sizeof(double));
  double *rotation = malloc(k*N * sizeof(double));
  
  for (int i=0; i<M*N; i++) x[i] = (double)i/M*N;
  
  BENCHMARK(pcapack_prcomp_svd(false, false, false, M, N, x, sdev, rotation), REPS);
  BENCHMARK(pcapack_prcomp_svd(true, false, false, M, N, x, sdev, rotation), REPS);
  BENCHMARK(pcapack_prcomp_svd(false, true, false, M, N, x, sdev, rotation), REPS);
  BENCHMARK(pcapack_prcomp_svd(true, true, false, M, N, x, sdev, rotation), REPS);
  
  free(x);
  free(sdev);
  free(rotation);
  return 0;
}
