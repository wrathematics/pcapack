#ifndef __PCAPACK_OMP_H__
#define __PCAPACK_OMP_H__


#define OMP_MIN_SIZE 2500


#ifdef _OPENMP
#include <omp.h>
#if _OPENMP >= 201307
#define _OPENMP_VER_4
#elif _OPENMP >= 200805
#define _OPENMP_VER_3
#endif
#endif


// Insert SIMD pragma if supported
#define SAFE_SIMD \
  #if _OPENMP_VER_4 \
    #pragma omp for simd \
  #endif

#endif
