! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


module variances
  implicit none
  
  
  contains
  
  subroutine variance(n, x, mn, var)
    ! in/out
    integer :: n
    double precision :: x( * ), mn, var
    ! local
    integer :: i
    double precision :: dt
    
    mn = 0
    var = 0
    
    do i = 1, n, 1
      dt = x(i) - mn
      mn = mn + dt/i
      var = var + dt*(x(i) - mn)
    end do
    
    var = var/(n-1)
    
    return
  end subroutine
  
  
  
  ! standard deviation
  subroutine standard_deviation(n, x, mn, sd)
    ! in/out
    integer :: n
    double precision :: x( * ), mn, sd
    ! subroutines
    external            dvar
    
    call dvar(n, x, mn, sd)
    sd = sqrt(sd)
    
    return
  end subroutine
  
  
  
  ! rowvars
  subroutine row_variances(m, n, x, ret)
    ! in/out
    integer :: m, n
    double precision :: x(m, n), ret( * )
    ! local
    integer :: i
    double precision :: mn, var
    ! subroutines
    external            dvar
    
    do i = 1, m, 1
      call dvar(n, x(i, :), mn, var)
      ret(i) = var
    end do
    
    return
  end subroutine
  
  
  
  ! colvars
  subroutine col_variances(m, n, x, ret)
    ! in/out
    integer :: m, n
    double precision :: x(m, n), ret( * )
    ! local
    integer :: i
    double precision :: mn, var
    ! subroutines
    external            dvar
    
    do i = 1, n, 1
      call dvar(n, x(:, i), mn, var)
      ret(i) = var
    end do
    
    return
  end subroutine
  
  
  
  ! rowvars
  subroutine row_standard_deviations(m, n, x, ret)
    ! in/out
    integer :: m, n
    double precision :: x(m, n), ret( * )
    ! local
    integer :: i
    double precision :: mn, var
    ! subroutines
    external            dvar
    
    do i = 1, m, 1
      call dvar(n, x(i, :), mn, var)
      ret(i) = sqrt(var)
    end do
    
    return
  end subroutine
  
  
  
  ! colvars
  subroutine col_standard_deviations(m, n, x, ret)
    ! in/out
    integer :: m, n
    double precision :: x(m, n), ret( * )
    ! local
    integer :: i
    double precision :: mn, var
    ! subroutines
    external            dvar
    
    do i = 1, n, 1
      call dvar(n, x(:, i), mn, var)
      ret(i) = sqrt(var)
    end do
    
    return
  end subroutine
  
  
end module
