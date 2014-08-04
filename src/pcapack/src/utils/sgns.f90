! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


! returns sign of x
module sgns
  implicit none
  
  
  interface sgn
    module procedure isgn
    module procedure ssgn
    module procedure dsgn
  end interface
  
  
  public :: sgn
  private :: isgn, ssgn, dsgn
  
  contains
  
  
  function isgn(x) result(xsign)
    integer :: xsgn
    integer, intent(in) :: x
    
    if (x > 0) then
      xsign = 1
    else if (x < 0) then
      xsign = -1
    else
      xsign = 0
    end if
    
    return
  end
  
  
  
  function ssgn(x) result(xsign)
    real :: xsign
    real, intent(in) :: x
    
    if (x > 0.0) then
      xsign = 1.0
    else if (x < 0) then
      xsign = -1.0
    else
      xsign = 0.0
    end if
    
    return
  end
  
  
  
  function dsgn(x) result(xsign)
    double precision :: xsign
    double precision, intent(in) :: x
    
    if (x > 0.0d0) then
      xsign = 1.0d0
    else if (x < 0) then
      xsign = -1.0d0
    else
      xsign = 0.0d0
    end if
    
    return
  end
  
end module


