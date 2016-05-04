#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "../misc.h"
#include "lapack.h"


// m == n

/**
 * @file
 * @brief Eigenvalue decomposition.
 *
 * @details
 * Eigenvalue decomposition.  Uses LAPACK's dsyevd. 
 *
 * @param inplace
 * Input.  Should the svd work on allocated x matrix, or a local copy?
 * @param only_values
 * Input.  Should only the values (not vectors) be computed?
 * @param symmetric
 * Input.  Is x symmetric?  In this case, only the upper triangle is used.
 * @maram n
 * Input.  Problem size (number of rows and columns of x).
 * @param x
 * In/Output data matrix.  If inplace==true then on successful function
 * return, x contains the eigenvectors.  Otherwise it's not modified (and
 * the vectors parameter contains the vectors).
 * @param values
 * Output.  The eigenvalues.
 * @param vectors
 * Output.  The eigenvectors (if inplace==false).
 *
 * @return
 * The return value indicates that status of the function.  Non-zero values
 * are errors.
 */
int pcapack_eig(const bool inplace, const bool only_values, const bool symmetric, const int n, double *restrict x, double *restrict values, double *restrict vectors)
{
  int info = 0;
  double *x_cp;
  // Passed to Fortran
  double worksize;
  int lwork, liwork;
  int *iwork;
  double *work;
  static int neg1 = -1;
  
  
  /*if (!inplace)*/
  /*{*/
    /*x_cp = malloc(n*n * sizeof(*x_cp));*/
    /*memcpy(x_cp, x, n*n*sizeof(double));*/
  /*}*/
  /*else*/
    /*x_cp = x;*/
  if (inplace)
    x_cp = x;
  else
  {
    memcpy(vectors, x, n*n*sizeof(double));
    x_cp = vectors;
  }
  
  
  if (likely(symmetric))
  {
    char uplo = 'u';
    char jobz;
    
    if (only_values)
      jobz = 'n';
    else
      jobz = 'v';
    
    dsyevd_(&jobz, &uplo, &n, x_cp, &n, values, &worksize, &neg1, &liwork, &neg1, &info);

    lwork = (int) worksize;
    work = malloc(lwork * sizeof(*work));
    iwork = malloc(liwork * sizeof(*iwork));
    
    dsyevd_(&jobz, &uplo, &n, x_cp, &n, values, work, &lwork, iwork, &liwork, &info);

    free(work);
    free(iwork);
  }
  #if 0 // TODO
  else
  {
    char jobvl = 'n', jobvr;
    
    if (only_values)
      jobvr = 'v';
    else
      jobvr = 'n';
    
    
    dgeev_(jobvl, jobvr, &n, x, &m, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info);
    
    
    dgeev_(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, &neg1, info);
  }
  #endif
  
  cleanup:
    if (!inplace) free(x_cp);

  return info;
}
/*dsyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info);*/


