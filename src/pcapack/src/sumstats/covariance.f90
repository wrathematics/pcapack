! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


module covariance
  use :: rank
  use :: sweeps
  use :: covariance_utils
  use :: variances
  use :: lapack, only : dlacpy
  implicit none
  
  
  integer(kind=C_int), parameter, public :: cor_pearson =  1
  integer(kind=C_int), parameter, public :: cor_spearman = 2
  integer(kind=C_int), parameter, public :: cor_kendall =  3 ! FIXME broken
  
  contains
  
  
  ! var-cov matrix of a matrix
  subroutine cov(m, n, x, ret)
    ! in/out
    integer :: m, n
    double precision :: x(m, n), ret(n, n)
    ! local
    integer :: i, j
    double precision :: alpha
    double precision, allocatable :: cpx(:, :)
    
    
    print *, m, n
    allocate(cpx(m, n))
!    call dlacpy('a', m, n, x, m, cpx, m)
    cpx = x
    
    ! center
    call center(m, n, cpx)
    
    ! compute variance-covariance matrix
    alpha = 1.0d0/(dble(m-1))
    
    call crossprod('n', m, n, alpha, cpx, ret)
    
    deallocate(cpx)
    
    return
  end
  
  
  
  ! correlation matrix
  subroutine cor(method, m, n, x, ret)
    ! in/out
    integer :: method
    integer :: m, n
    double precision :: x(m, n), ret(n, n)
    ! local
    integer :: i, j
    double precision, allocatable :: rnks(:, :)
    
    
    if (method == cor_pearson) then
!      call variance(m, n, x, ret)
!      call scaler(n, n, ret)
    else if (method == cor_spearman) then
      allocate(rnks(m, n))
      
      call colrank(rank_min, m, n, x, ret)
      call cov(m, n, rnks, ret)
      
      deallocate(rnks)
!      else if (method == rank_kendall) then
!        
    end if
    
    return
  end
  
end module

