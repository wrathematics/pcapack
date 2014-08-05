! this source code form is subject to the terms of the mozilla public
! license, v. 2.0. if a copy of the mpl was not distributed with this
! file, you can obtain one at http://mozilla.org/mpl/2.0/.

! copyright 2013, schmidt

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
!  scale_      (input) character*1
!           flag for whether or not the data should first be scaled
!           by the column standard deviations.  if scale_='y' then the
!           data is scaled, and otherwise not.
!
!  info     (output) integer
!           lapack-style error return.  

module pca
  use :: lapack, only : dgemm, dscal, dlacpy, dgesdd
  use, intrinsic :: iso_c_binding
  use :: sweeps
  implicit none
  
  
  logical(kind=c_bool), parameter, public :: true = .true.
  logical(kind=c_bool), parameter, public :: false = .false.
  
  interface operator (==)
    module procedure cbool_equality
  end interface
  
  contains
  
  
  function cbool_equality(v1,v2) result (v3)
    logical(kind=c_bool), intent(in) :: v1, v2
    logical :: v3
    v3 = v1 .eqv. v2
    return
  end function
  
  
  subroutine prcomp_svd_work(m, n, k, x, sdev, trot, retrot, center_, scale_, info)
    ! in/out
    logical(kind=C_bool), intent(in) :: retrot, center_, scale_
    integer, intent(in) :: m, n, k
    integer, intent(out) :: info
    double precision, intent(inout) :: x(m, n), sdev(k), trot(k, n)
    ! local
    integer :: lwork
    double precision :: tmp(1)
    double precision, allocatable :: work(:), u(:)
    integer,          allocatable :: iwork(:)
    
    
    ! allocations
    allocate(iwork(8*k))
    lwork = -1
    
    call dgesdd('o', m, n, x, m, sdev, u, m, trot, k, tmp, lwork, iwork, info)
   
    lwork = int(tmp(1))
    
    allocate(work(lwork))
    
    if (center_ == true .and. scale_ == true) then
      call center_scale(m, n, x)
    else if (center_ == true) then
      call center(m, n, x)
    else if (scale_ == true) then
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
  
  
  
  subroutine prcomp_svd(m, n, k, x, sdev, trot, retrot, center_, scale_, info) &
    bind(C, name='prcomp_svd_')
    ! in/out
    logical(kind=c_bool), intent(in) :: retrot, center_, scale_
    integer :: m, n, k, info
    double precision :: x(m, n), sdev(k), trot(k, n)
    ! local
    double precision, allocatable :: cpx(:,:)
    
    
    ! allocations
    allocate(cpx(m, n))
    
    ! compute pc's
    if (retrot .eqv. true) then ! return rotated x
      call dlacpy('a', m, n, x, m, cpx, m)
      
      call prcomp_svd_work(m, n, k, x, sdev, trot, retrot, center_, scale_, info)
      
      ! x = x %*% t(rot)            c := alpha*op(a)*op(b) + beta*c
      call dgemm('n', 't', m, n, k, 1.0d0, cpx, m, trot, k, 0.0d0, x, m)
      
      deallocate(cpx)
    else
      call prcomp_svd(m, n, k, x, sdev, trot, retrot, center_, scale_, info)
    end if
    
    
    return
  end subroutine
  
  
  
! principal components via eigenvalue decomposition of covariance matrix
!  subroutine prcomp_eig(m, n, k, x, sdev, trot, retrot, center, scale_, info)
!    
!    return
!  end subroutine
  
end module
