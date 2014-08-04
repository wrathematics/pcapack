! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


module distances
  use :: lapack
  implicit none
  
  integer, public, parameter :: dist_euclidean = 1
  integer, public, parameter :: dist_supremum  = 2
  integer, public, parameter :: dist_infimum   = 3
  integer, public, parameter :: dist_manhattan = 4
  integer, public, parameter :: dist_rcanberra = 5
  integer, public, parameter :: dist_minkowski = 6
  integer, public, parameter :: dist_canberra  = 7
  !  integer, public, parameter :: dist_binary    = 8
  
  
  contains
  
  
  function distance(method, n, x, y, p) result(dist)
    double precision :: dist
    ! in/out
    integer, intent(in) :: method, n, p
    double precision, intent(in) :: x(n), y(n)
    ! local
    integer :: i, allocerr
    integer :: ux, uy
    double precision, allocatable :: work(:)
    
    
    allocerr = 0
    allocate(work(n), stat=allocerr)
    if (allocerr /= 0) stop "out of memory"
    
    ! work = x - y
    call dcopy(n, x, 1, work, 1)
    
    call daxpy(n, -1.0d0, y, 1, work, 1)
    
    if (method == dist_euclidean) then
      dist = dnrm2(n, work, 1)
    else if (method == dist_supremum) then
      dist = work(idamax(n, work, 1))
    else if (method == dist_infimum) then
      dist = work(idamin(n, work, 1))
      return
    else if (method == dist_manhattan) then
      dist = dasum(n, work, 1)
    else if (method == dist_rcanberra) then
      dist = 0.0d0
      do i = 1, n
        dist = dist + abs(work(i)) / abs(x(i-1) + y(i-1))
      end do
    else if (method == dist_canberra) then
      dist = 0.0d0
      do i = 1, n
        dist = dist + abs(work(i)) / (abs(x(i-1)) + abs(y(i-1)))
      end do
    else if (method == dist_minkowski) then
      dist = dnrm3(n, work, 1, p)
!    else if (method == dist_binary) then
!      dist = 0.0d0
!      do i = 1, n
!        ux = transfer(x(i), 1)
!        uy = transfer(y(i), 1)
!        
!        
!      end do
    end if
    
    deallocate(work)
    
    return
  end function
  
  
  
  ! mahalanobis distance d_m(x,y)
  function distance_mahalanobis(n, x, y, covinv) result(mahal)
    double precision :: mahal
    ! in/out
    integer, intent(in) :: n
    double precision, intent(in) :: x(n), y(n), covinv(*)
    ! local
    integer :: i, allocerr
    double precision, allocatable :: work1(:)
    double precision, allocatable :: work2(:)
    
    
    allocerr = 0
    allocate(work1(n), stat=allocerr)
    if (allocerr.ne.0) stop "out of memory"
    allocate(work2(n), stat=allocerr)
    if (allocerr.ne.0) stop "out of memory"
    
    ! work1 = x - y
    call dcopy(n, x, 1, work1, 1)
    
!    call daxpy(n, 1.0d0, work1, 1, y, 1)
    
    ! mahal = work1^t * covinv * work2
    call dgemv('n', n, n, 1.0d0, covinv, n, work1, 1, 0.0d0, work2, 1)
    
    mahal = sqrt(ddot(n, work1, 1, work2, 1))
    
    deallocate(work1)
    deallocate(work2)
    
    return
  end function
  
  
  
!  subroutine distance_matrix(method, m, n, x, ret, p)
!    integer, intent(in) :: method, n, p
!    double precision, intent(in) :: x(m, n), ret(*)
!    integer :: retlen
!    
!    
!    retlen = n*(n-1)/2
!    
!    
!    
!  end subroutine
end module

