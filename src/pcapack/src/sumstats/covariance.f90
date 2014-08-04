! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


module covariance
  use :: covariance_utils
  implicit none
  
  
  contains
  
  
  ! var-cov matrix of a matrix
  pure subroutine cov(m, n, x, cov, info)
    ! in/out
    integer             m, n, info
    double precision    x(m, n), cov(n, n)
    ! local
    integer             i, j
    double precision    alpha
    double precision, allocatable :: cpx(:, :)
    ! external
    external           dlacpy, dcrossprod, dcntr1, dinvip
    
    
    allocate(cpx(m, n))
    call dlacpy('a', m, n, x, m, cpx, m)
    
    ! center
    call center(m, n, cpx)
    
    ! compute variance-covariance matrix
    alpha = 1/(dble(m)-1)
    
    call crossprod('n', m, n, alpha, cpx, cov)
    
    deallocate(cpx)
    
    return
  end
  
  
  
  ! correlation matrix
  pure subroutine cor(method, m, n, x, cor)
    ! in/out
    character*1           method
    integer             m, n
    double precision    x(m, n), cor(n, n)
    ! local
    integer             i, j
    double precision, allocatable :: rnks(:, :)
    ! external
    external           dcov
    
    
    if (method == 'p' .or. method == 'P') then
      call dvar(m, n, x, cor)
      call dscl(n, n, cor) ! fixme
    else if (method == 's' .or. method == 'S') then
      allocate(rnks(m, n))
      
      call dclrank(m, n, x, rnks, 1)
      call dcov(m, n, rnks, cor)
      
      deallocate(rnks)
!      else if (method == 'k') then
!        
    end if
    
    return
  end
  
  
end module
