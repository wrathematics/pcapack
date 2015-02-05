// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Copyright 2015, Schmidt

#include <math.h>

#define MIN(x,y) x<y?x:y
#define MAX(x,y) x<y?y:x

// TODO: xpose, center_scale, center, scale

int pcapack_prcomp_svd(int m, int n, double *x, double *sdev, double *rotation, bool retrot, bool centerx, bool scalex)
{
  char trans = 'n';
  int info;
  int minmn = MIN(m, n);
  int intone = 1;
  double *x_cp;
  double *u;
  double one = 1., zero = 0.;
  double tmp;
  
  
  x_cp = malloc(m*n * sizeof(*x));
  memcpy(x_cp, x, m*n*sizeof(*x));
  
  u = malloc(m*minmn * sizeof(*x));
  
  if (centerx && scalex)
    center_scale(m, n, x_cp);
  else if (centerx)
    center(m, n, x_cp);
  else if (scalex)
    scale(m, n, x_cp);
  
  info = pcapack_svd(n, n, m, n, x_cp, sdev, u, rotation)
  
  // TODO check info
  
  xpose(minmn, n, rotation);
  
  if (retrot)
    dgemm_(&trans, &trans, &m, &minmn, &n, &one, x, &m, rotation, &n, &zero, x, &m);
  
  tmp = 1. / MAX(1., sqrt((double) m-1));
  dscal_(&minmn, &tmp, sdev, &intone);
  
  free(x_cp);
  free(u);
  
  return info;
}
    ! normalize singular values
    tmp = 1.0d0 / max(1.0d0, dsqrt(dble(m-1)))
    call dscal(minmn, tmp, sdev, 1)
    
    
    if (allocated(u)) deallocate(u)
    if (allocated(cpx)) deallocate(cpx)
    
    return
  end subroutine
  
  


/*
  subroutine prcomp_eig(m, n, x, sdev, rotation, retrot, info) &
  bind(C, name='prcomp_eigcov_')
    ! in/out
    logical(kind=c_bool), intent(in) :: retrot
    integer, intent(in) :: m, n
    integer, intent(out) :: info
    double precision, intent(inout) :: x(m, n)
    double precision, intent(out) :: sdev(*), rotation(*)
    ! local
    integer :: i
    integer :: iworksize(1), lwork, liwork
    double precision :: tmp, worksize(1)
    double precision, allocatable :: covmat(:,:), work(:), iwork(:)
    
    
    allocate(covmat(n, n))
    
    call cov(m, n, x, covmat)
    
    tmp = 1.0d0 - 1.0d0/dble(m)
    call dscal(n*n, tmp, covmat, 1)
    
    ! compute eigen
    call dsyevd('V', 'U', n, covmat, n, sdev, worksize, -1, iworksize, -1, info)
    lwork = int(worksize(1))
    liwork = iworksize(1)
    allocate(work(lwork))
    allocate(iwork(liwork))
    call dsyevd('V', 'U', n, covmat, n, sdev, work, lwork, iwork, liwork, info)
    
    ! sdev = rev(sqrt(sdev))
    do i = 1, n/2
      tmp = sdev(i)
      sdev(i) = dsqrt(sdev(n-i+1))
      sdev(n-i+1) = dsqrt(tmp)
    end do
    
    if (1 == int(mod(n, 2))) sdev(n/2+1) = sqrt(sdev(n/2+1))
    
    
    call xpose(n, n, rotation)
    
    if (retrot == true) then
      call dgemm('n', 'n', m, n, n, 1.0d0, x, m, rotation, n, 0.0d0, x, m)
    end if
    
    
    if (allocated(covmat)) deallocate(covmat)
    if (allocated(work)) deallocate(work)
    if (allocated(iwork)) deallocate(iwork)
    
    return
  end subroutine
*/
