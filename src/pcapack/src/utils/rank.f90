! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


module rank
  implicit none
  
  
  interface
    pure subroutine genrank(x, rank, length, method) bind(C)
      use, intrinsic :: iso_c_binding
      real(kind=c_double), intent(in) :: x
      real(kind=c_double), intent(out) :: rank
      integer, intent(in), value :: length, method
    end subroutine
  end interface
  
  
  contains
  
  
  pure subroutine colrank(m, n, x, ret, method)
    ! In/Out
    integer, intent(in) :: m, n, method
    double precision, intent(in) :: x(m, n)
    double precision, intent(out) :: ret(m, n)
    integer :: i, j
    
    do i = 1, n, 1
      call genrank(x(1, i), ret(1, i), m, method)
    end do
    
    return
  end subroutine
  
end module

