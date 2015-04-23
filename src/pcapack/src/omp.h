#ifndef __PCAPACK_OMP_H__
#define __PCAPACK_OMP_H__


#define OMP_MIN_SIZE 2500


#ifdef _OPENMP
#include <omp.h>
#if _OPENMP >= 201307
#define OMP_VER_4
#elif _OPENMP >= 200805
#define OMP_VER_3
#endif
#endif


// Insert SIMD pragma if supported
#ifdef OMP_VER_4
#define SAFE_SIMD _Pragma("omp simd")
#define SAFE_FOR_SIMD _Pragma("omp for simd")
#else
#define SAFE_SIMD 
#define SAFE_FOR_SIMD
#endif


// Inform compiler of byte alignment if GCC >= 4.7
#ifdef __GNUC__
#if __GNUC__ >= 4 && __GNUC_MINOR__ >= 7
#define ALIGNMENT(x,bytes) x=__builtin_assume_aligned(x,bytes)
#else
#define ALIGNMENT(x,bytes)
#endif
#endif



#endif
