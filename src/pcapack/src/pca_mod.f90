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

module pca_mod
  use :: lapack
  use, intrinsic :: iso_c_binding
  use :: sweeps
  use :: cbool
  use :: svd_mod
  use :: transposition
  use :: covariance
  implicit none
  
  
  contains
  
  subroutine prcomp_svd(m, n, x, sdev, rotation, retrot, centerx, scalex, info) &
    bind(C, name='prcomp_svd_')
    ! in/out
    logical(kind=c_bool), intent(in) :: retrot, centerx, scalex
    integer, intent(in) :: m, n
    integer, intent(out) :: info
    double precision, intent(inout) :: x(m, n)
    double precision, intent(out) :: sdev(*), rotation(*)
    ! local
    integer :: minmn
    double precision :: tmp
    double precision, allocatable :: u(:,:), cpx(:,:)
    
    
    minmn = min(m, n)
    
    ! allocations
    allocate(cpx(m, n))
    
    call dlacpy('a', m, n, x, m, cpx, m)
    
    allocate(u(m, minmn))
    
    
    ! center/scale
    if (centerx == true .and. scalex == true) then
      call center_scale(m, n, cpx)
    else if (centerx == true) then
      call center(m, n, cpx)
    else if (scalex == true) then
      call scaler(m, n, cpx)
    end if
    
    
    ! compute svd
    call svd(n, n, m, n, cpx, sdev, u, rotation, info)
    
    call xpose(minmn, n, rotation)
    
    if (retrot == true) then
      call dgemm('n', 'n', m, minmn, n, 1.0d0, x, m, rotation, n, 0.0d0, x, m)
    end if
    
    
    ! normalize singular values
    tmp = 1.0d0 / max(1.0d0, dsqrt(dble(m-1)))
    call dscal(minmn, tmp, sdev, 1)
    
    
    if (allocated(u)) deallocate(u)
    if (allocated(cpx)) deallocate(cpx)
    
    return
  end subroutine
  
  
  
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
  
end module
