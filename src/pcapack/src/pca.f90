! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


! purpose
! =======
! this routine computes the principal component, 
!
!
! notes
! =====
!
! this method computes the principal components decomposition of
! a given input matrix x.
!
!
! arguments
! =========
!
!  m, n     (input) integer
!           the number of rows and columns of the input matrix a.
!
!  a        (input) double precision
!           input data matrix of size mxn.
!
!  k        (input) integer
!           the desired number of singular vectors.
!
!  sdev     (output) double precision
!           the standard deviations of the principal components; 
!           equivalently, the singular values.
!
!  trot     (output) double precision
!           
!
!  retrot   (input) character*1
!           flag for whether or not the rotated values should be 
!           returned, i.e., if trot=x*v^t should be computed, where 
!           v is the matrix of right singular vectors.  if 
!           retrot='y' then trot is computed, and otherwise ignored.
!
!  center     (input) character*1
!           flag for whether or not the data should first be centered.
!           if center='y' then the data is first centered, and otherwise
!           not.
!
!  scalex      (input) character*1
!           flag for whether or not the data should first be scaled
!           by the column standard deviations.  if scalex='y' then the
!           data is scaled, and otherwise not.
!
!  info     (output) integer
!           lapack-style error return.  

module pca
  use :: lapack
  use, intrinsic :: iso_c_binding
  use :: sweeps
  use :: cbool
  implicit none
  
  
  contains
  
  subroutine prcomp_svd_work(m, n, k, x, sdev, trot, retrot, centerx, scalex, info)
    ! in/out
    logical(kind=C_bool), intent(in) :: retrot, centerx, scalex
    integer, intent(in) :: m, n, k
    integer, intent(out) :: info
    double precision, intent(inout) :: x(m, n), sdev(k), trot(k, n)
    ! local
    integer :: lwork
    double precision :: tmp(1)
    double precision, allocatable :: work(:), u(:)
    integer, allocatable :: iwork(:)
    
    
    ! allocations
    allocate(iwork(8*k))
    lwork = -1
    
    call dgesdd('o', m, n, x, m, sdev, u, m, trot, k, tmp, lwork, iwork, info)
   
    lwork = int(tmp(1))
    
    allocate(work(lwork))
    
    if (centerx == true .and. scalex == true) then
      call center_scale(m, n, x)
    else if (centerx == true) then
      call center(m, n, x)
    else if (scalex == true) then
      call scaler(m, n, x)
    end if
      
    ! compute svd
    if (m < n) allocate(u(m*m))
    call dgesdd('o', m, n, x, m, sdev, u, m, trot, k, work, lwork, iwork, info)
    if (m < n) deallocate(u)
    
    deallocate(iwork)
    deallocate(work)
    
    ! normalize singular values
    tmp(1) = 1.0d0 / max(1.0d0, sqrt(dble(m-1)))
    call dscal(k, tmp(1), sdev, 1)
    
    return
  end subroutine
  
  
  
  subroutine prcomp_svd(m, n, k, x, sdev, trot, retrot, centerx, scalex, info) &
    bind(C, name='prcomp_svd_')
    ! in/out
    logical(kind=c_bool), intent(in) :: retrot, centerx, scalex
    integer :: m, n, k, info
    double precision :: x(m, n), sdev(k), trot(k, n)
    ! local
    double precision, allocatable :: cpx(:,:)
    
    
    ! allocations
    allocate(cpx(m, n))
    
    call dlacpy('a', m, n, x, m, cpx, m)
    ! compute pc's
    call prcomp_svd_work(m, n, k, cpx, sdev, trot, retrot, centerx, scalex, info)
      
!    if (retrot == true) then ! return rotated x
!      ! x = x %*% t(rot)            c := alpha*op(a)*op(b) + beta*c
!      call dgemm('n', 't', m, n, k, 1.0d0, cpx, m, trot, k, 0.0d0, cpx, m)
!    end if
    
    deallocate(cpx)
    
    return
  end subroutine
  
  
  
! principal components via eigenvalue decomposition of covariance matrix
!  subroutine prcomp_eig(m, n, k, x, sdev, trot, retrot, center, scalex, info)
!    
!    return
!  end subroutine
  
end module
