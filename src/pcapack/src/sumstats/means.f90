! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


module means
  implicit none
  
  contains
  
  
  subroutine rowsums(m, n, x, ret)
    ! in/out
    integer :: m, n
    double precision :: x(m, n), ret(m)
    integer :: i
    
    
    do i = 1, m
      ret(i) =  sum(x(i, :))
    end do
    
    return
  end subroutine
  
  
  
  subroutine colsums(m, n, x, ret)
    ! in/out
    integer :: m, n
    double precision :: x(m, n), ret(n)
    ! local
    integer :: j
    
    do j = 1, n
      ret(j) = sum(x(:, j))
    end do
    
    return
  end subroutine
  
  
  
  double precision function mean(n, x)
    ! in
    integer :: n
    double precision :: x(n)
    ! local
    integer :: i
    
    mean = 0
    
    do i = 1, n
      mean = mean + x(i)/n
    end do
    
    return
  end function
  
  
  
  subroutine row_means(m, n, x, ret)
    implicit none
    ! in/out
    integer :: m, n
    double precision :: x(m, n), ret( * )
    ! local
    integer :: i, j
    
    do j = 1, n
      do i = 1, m
        ret(i) = ret(i) + x(i, j) / n
      end do
    end do
    
    return
  end subroutine
  
  
  
  subroutine col_means(m, n, x, ret)
    implicit none
    ! in/out
    integer :: m, n
    double precision :: x(m, n), ret( * )
    ! local
    integer :: i, j
    
    do j = 1, n
      do i = 1, m
        ret(j) = ret(j) + x(i, j) / m
      end do
    end do
    
    return
  end subroutine
  
end module
