! this source code form is subject to the terms of the mozilla public
! license, v. 2.0. if a copy of the mpl was not distributed with this
! file, you can obtain one at http://mozilla.org/mpl/2.0/.

! copyright 2013, schmidt

! distance d_{method}(x,y)
! methods:  'e' - euclidean
!           's' - supremum (infinity)
!           'i' - infimum (-infinity)
!           'h' - manhattan
!           'c' - canberra
!           'k' - minkowski

module distances
  implicit none
  
  integer, public, parameter :: dist_euclidean = 1
  integer, public, parameter :: dist_supremum  = 2
  integer, public, parameter :: dist_infimum   = 3
  integer, public, parameter :: dist_manhattan = 4
  integer, public, parameter :: dist_canberra  = 5
  integer, public, parameter :: dist_minkowski = 6
  
  
  contains
  
  
  function distance(method, n, x, y, p) 
    implicit none
    ! in/out
    integer, intent(in) :: method, n, p
    double precision, intent(in) :: x(n), y(n)
    ! local
    integer :: i, allocerr
    double precision, allocatable :: work(:)
    ! functions
    external            dcopy, daxpy
    integer :: idamax, idamin
    double precision :: dasum, dnrm2, dnrm3
    
    
    allocerr = 0
    allocate(work(n), stat=allocerr)
    if (allocerr.ne.0) stop "out of memory"
    
    ! work = x - y
    call dcopy(n, x, 1, work, 1)
    
    call daxpy(n, 1.0d0, work, 1, y, incy)
    
    if (method == dist_euclidean) then
      ddist = dnrm2(n, work, 1)
    else if (method == dist_supremum) then
      ddist = x(idamax(n, x, 1))
    else if (method == dist_infimum) then
      ddist = x(idamin(n, x, 1))
    else if (method == dist_manhattan) then
      ddist = dasum(n, work, 1)
    else if (method == dist_canberra) then
      ddist = 0.0d0
      do i = 1, n
        ddist = ddist + abs(work(i)) / (abs(x((i-1) + 1)) + abs(y((i-1) + 1)))
      end do
    else if (method == dist_minkowski) then
      ddist = dnrm3(n, work, 1, p)
    end if
    
    deallocate(work)
    
    return
  end function
  
  
  
  ! mahalanobis distance d_m(x,y)
  function distance_mahalanobis(n, x, y, covinv) 
  result(dmahal)
  
    implicit none
    ! in/out
    integer, intent(in) :: n
    double precision, intent(in) :: x(n), y(n), covinv(*)
    ! local
    integer :: i, allocerr
    double precision, allocatable :: work1(:)
    double precision, allocatable :: work2(:)
    ! external
    external            dcopy, daxpy, dgemv
    double precision :: ddot
    
    
    allocerr = 0
    allocate(work1(n), stat=allocerr)
    if (allocerr.ne.0) stop "out of memory"
    allocate(work2(n), stat=allocerr)
    if (allocerr.ne.0) stop "out of memory"
    
    ! work1 = x - y
    call dcopy(n, x, 1, work1, 1)
    
    call daxpy(n, 1.0d0, work1, 1, y, 1)
    
    ! dmahal = work1^t * covinv * work2
    call dgemv('n', n, n, 1.0d0, covinv, n, work1, 1, 0.0d0, work2, 1)
    
    dmahal = sqrt(ddot(n, work1, 1, work2, 1))
    
    deallocate(work1)
    deallocate(work2)
    
    return
  end function

end module

