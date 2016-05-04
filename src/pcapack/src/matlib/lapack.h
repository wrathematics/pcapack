/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2015, Schmidt

#ifndef __PCAPACK_LAPACK_H__
#define __PCAPACK_LAPACK_H__


// BLAS
void dscal_(int *n, double *a, double *x, int *incx);
void daxpy_(int *n, double *a, double *x, int *incx, double *y, int *incy);

void dgemv_(const char *trans, const int *m, const int *n, const double *restrict alpha, 
            const double *restrict a, const int *lda, const double *restrict x, const int *incx, 
            const double *restrict beta, double *restrict y, const int *incy);
void dger_(const int *m, const int *n, const double *restrict alpha, const double *restrict x, 
           const int *incx, const double *restrict y, const int *incy, double *restrict a, 
           const int *lda);

void dsyrk_(const char *uplo, const char *trans, const int *n, const int *k, 
            const double *restrict alpha, const double *restrict a, const int *lda, 
            const double *restrict beta, double *restrict c, const int *ldc);
void dgemm_(const char *transa, const char *transb, const int *m, const int *n, 
            const int *k, const double *restrict alpha, const double *restrict a, 
            const int *lda, const double *restrict b, const int *ldb, 
            const double *beta, double *restrict c, const int *ldc);



// LAPCK
void dlassq_(int *n, double *x, int *incx, double *scale, double *sumsq);

void dgeqrf_(const int *m, const int *n, double *restrict a, const int *lda, 
             double *restrict tau, double *restrict work, const int *restrict lwork, 
             int *info);
void dorgqr_(const int *m, const int *n, const int *k, double *restrict a, 
             const int *lda, double *restrict tau, double *restrict work, const int *lwork, int *info);

void dgetrf_(const int *m, const int *n, double *restrict a, const int *lda, int *ipiv, int *info);
void dgetri_(const int *n, double *restrict a, const int *lda, int *ipiv, double *restrict work, int *lwork, int *info);

void dgesdd_(const char *jobz, const int *m, const int *n, double *a, 
             const int *lda, double *restrict s, double *restrict u, const int *restrict ldu, double *restrict vt, 
             const int *ldvt, double *restrict work, const int *lwork, int *iwork, int *info);
void dsyevd_(const char *jobz, const char *uplo, const int *n, double *restrict a, 
             const int *lda, double *restrict w, double *restrict work, const int *lwork, int *restrict iwork, 
             const int *liwork, int *info);
void dsyevr_(char *jobz, char *range, char *uplo, int *n, double *a, int *lda, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w, double *z, int *ldz, int *isuppz, double *work, int *lwork, int *iwork, int *liwork, int *info);
void dgeev_(char *jobvl, char *jobvr, int *n, double *a, int *lda, double *wr, double *wi, double *vl, int *ldvl, double *vr, int *ldvr, double *work, int *lwork, int *info);



// Custom
double dnrm2(int n, double *x, int incx);


#endif
