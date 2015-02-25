#ifndef __PCAPACK_BENCHMARK_H__
#define __PCAPACK_BENCHMARK_H__


#include <stdio.h>
#include <time.h>
#include <sys/time.h>

static double systime_wall()
{
  struct timeval t;
  gettimeofday(&t, NULL);
  return (double)t.tv_sec + ((double) t.tv_usec) / 1000000.;
}

static double systime_cpu()
{
    return (double) clock() / CLOCKS_PER_SEC;
}

double __benchmark_wall_seconds, __benchmark_cpu_seconds;
double __benchmark_wall_start, __benchmark_wall_end;
double __benchmark_cpu_start, __benchmark_cpu_end;

#define MAXPRINTLEN 20
static void print_func_name(char *call)
{
  int i = 0;
  putchar('\n');
  while (call[i] != '(' && call[i] != '\0' && i<MAXPRINTLEN)
    putchar(call[i++]);
  
  printf(":\n");
}

#define TIMEEXPR(foo) \
  __benchmark_cpu_start = systime_cpu(); \
  __benchmark_wall_start = systime_wall(); \
  foo; \
  __benchmark_cpu_end = systime_cpu(); \
  __benchmark_wall_end = systime_wall(); \
  __benchmark_cpu_seconds = __benchmark_cpu_end - __benchmark_cpu_start; \
  __benchmark_wall_seconds = __benchmark_wall_end - __benchmark_wall_start;

#define SHOWTIMES(foo) \
  print_func_name(#foo); \
  printf("\tcpu\twall\n"); \
  printf("Total:\t%.3f\t%.3f\n", __benchmark_cpu_seconds, __benchmark_wall_seconds);

#define SYSTIME(foo) \
  TIMEEXPR(foo); \
  SHOWTIMES(foo);

#define REPEAT(foo, n) \
  for (volatile int __benchint_i=0; __benchint_i<n; __benchint_i++) foo;

#define BENCHMARK(foo, n) \
  TIMEEXPR(REPEAT(foo, n)); \
  SHOWTIMES(foo); \
  printf("Avg:  \t%.3f\t%.3f\n", __benchmark_cpu_seconds/n, __benchmark_wall_seconds/n);


#endif
