! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt

! Inspired by http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/


module float_equals
  use :: sgns
  use :: lapack, only : dlamch
  implicit none
  
  
  contains
  
  function fequals(x, y) result(ret)
    logical :: ret
    double precision, intent(in) :: x, y
    double precision :: tol
    integer :: ulpsdiff, ux, uy
    integer :: ulpstol = 2
    intrinsic :: dsqrt, dabs, iabs
    
    
    tol = dsqrt(dlamch('e')) ! inherited from R
    
    ret = dabs(x-y) < tol
    if (ret) return
    
    
    if(sign(1.0d0, x) /= sign(1.0d0, y)) then
      ! Check for +0/-0 case
      if (x == y) then
        ret = .true.
      else
        ret = .false.
      end if
      
      return
    end if
    
    ! Difference in ULP's
    ux = transfer(x, 1)
    uy = transfer(y, 1)
    
    ulpsdiff = iabs(ux - uy)
    if (ulpsdiff <= ulpstol) then
      ret = .true.
      return
    end if
    
    
    ret = .false.
    
    return
  end function
  
end module

