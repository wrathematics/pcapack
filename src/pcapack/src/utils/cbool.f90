! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


module cbool
  use, intrinsic :: iso_c_binding
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
  
end module
