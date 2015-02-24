#ifndef __PCAPACK_BENCHMARK_H__
#define __PCAPACK_BENCHMARK_H__


#include <stdio.h>
#include <time.h>

double __benchmark_seconds;
clock_t __benchmark_start, __benchmark_end;

#define TIMEEXPR(foo) \
  __benchmark_start = clock(); \
  foo; \
  __benchmark_end = clock(); \
  __benchmark_seconds = (double)(__benchmark_end - __benchmark_start) / CLOCKS_PER_SEC

#define SYSTIME(foo) \
  TIMEEXPR(foo) \
  printf("%s:\t%d.%d seconds\n", #foo, __benchmark_msec/1000, __benchmark_msec%1000)

#define REPEAT(foo, n) \
  for (int __benchint_i=0; __benchint_i<n; __benchint_i++) foo;

#define BENCHMARK(foo, n) \
  TIMEEXPR(REPEAT(foo, n)); \
  printf("%s:\n", #foo); \
  printf("    %.3f seconds total\n", __benchmark_seconds); \
  printf("    %.3f seconds average (%d runs)\n", __benchmark_seconds/((double) n), n)


#endif
