! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


module sweeps
  implicit none
  
  
  contains
  
  
  integer function ind(i, j)
    integer i, j
    
    ind = mod(i, j)
    if (ind == 0) then
      ind = j
    end if
    
    return
  end function
  
  
  
  ! sweep array out of matrix
  ! code is a bit of a monstrosity, but not sure how to simplify it without hurting performance...
  subroutine sweep(m, n, x, vec, lvec, margin, fun)
    implicit none
    ! in/out
    integer :: m, n, lvec, margin
    double precision :: x(m,n), vec(lvec)
    character*1           fun
    ! local
    integer :: i, j, pos
    ! function
    integer :: ind
    
    ! special case 1
    if (margin == 1 .and. lvec == m) then
      if (fun == "+") then
        do i = 1, m
          x(i, :) = x(i, :) + vec(i)
        end do
      else if (fun == "-") then
        do i = 1, m
          x(i, :) = x(i, :) - vec(i)
        end do
      else if (fun == "*") then
        do i = 1, m
          x(i, :) = x(i, :) * vec(i)
        end do
      else if (fun == "/") then
        do i = 1, m
          x(i, :) = x(i, :) / vec(i)
        end do
      end if
    ! special case 2
    else if (margin == 2 .and. lvec == n) then
      if (fun == "+") then
        do j = 1, n
          x(:, j) = x(:, j) + vec(j)
        end do
      else if (fun == "-") then
        do j = 1, n
          x(:, j) = x(:, j) - vec(j)
        end do
      else if (fun == "*") then
        do j = 1, n
          x(:, j) = x(:, j) * vec(j)
        end do
      else if (fun == "/") then
        do j = 1, n
          x(:, j) = x(:, j) / vec(j)
        end do
      end if
    ! general case
    else
      ! addition
      if (fun == "+") then
        if (margin == 1) then
          pos = 1
          do j = 1, n
            do i = 1, m
              x(i, j) = x(i, j) + vec(pos)
              pos = ind(pos+1, lvec)
            end do
          end do
        else if (margin == 2) then
          do j = 1, n
            pos = ind(j, lvec)
            do i = 1, m
              x(i, j) = x(i, j) + vec(pos)
              pos = ind(pos+n, lvec)
            end do
            pos = ind(pos + 1, lvec)
          end do
        end if
      ! subtraction
      else if (fun == "-") then
        if (margin == 1) then
          pos = 1
          do j = 1, n
            do i = 1, m
              x(i, j) = x(i, j) - vec(pos)
              pos = ind(pos+1, lvec)
            end do
          end do
        else if (margin == 2) then
          do j = 1, n
            pos = ind(j, lvec)
            do i = 1, m
              x(i, j) = x(i, j) - vec(pos)
              pos = ind(pos+n, lvec)
            end do
            pos = ind(pos + 1, lvec)
          end do
        end if
      ! multiplication
      else if (fun == "*") then
        if (margin == 1) then
          pos = 1
          do j = 1, n
            do i = 1, m
              x(i, j) = x(i, j) * vec(pos)
              pos = ind(pos+1, lvec)
            end do
          end do
        else if (margin == 2) then
          do j = 1, n
            pos = ind(j, lvec)
            do i = 1, m
              x(i, j) = x(i, j) * vec(pos)
              pos = ind(pos+n, lvec)
            end do
            pos = ind(pos + 1, lvec)
          end do
        end if
      ! division
      else if (fun == "/") then
        if (margin == 1) then
          pos = 1
          do j = 1, n
            do i = 1, m
              x(i, j) = x(i, j) / vec(pos)
              pos = ind(pos+1, lvec)
            end do
          end do
        else if (margin == 2) then
          do j = 1, n
            pos = ind(j, lvec)
            do i = 1, m
              x(i, j) = x(i, j) / vec(pos)
              pos = ind(pos+n, lvec)
            end do
            pos = ind(pos + 1, lvec)
          end do
        end if
      end if
    end if
    
    return
  end subroutine
  
  
  
  ! centers matrix x; return overwrites x
  ! if centering and scaling is needed, call dcntrscl; it's slightly more efficient
  ! the subroutine passes through each column, determining the mean 
  subroutine center(m, n, x)
    implicit none
    ! in/out
    integer :: m, n
    double precision :: x(m, n)
    ! local
    integer :: j
    double precision :: mn
    
    do j = 1, n, 1
      mn = sum(x(:, j)/m)
      x(:, j) = x(:, j) - mn
    end do
    
    return
  end subroutine
  
  
  
  ! scales matrix x by root mean square; return overwrites x
  ! rms = sqrt(sum(x^2)/(n-1))
  ! if centering and scaling is needed, call dcntrscl; it's slightly more efficient
  subroutine scaler(m, n, x)
    implicit none
    ! in/out
    integer :: m, n
    double precision :: x(m, n)
    ! local
    integer :: i, j, p
    double precision :: scl, dt
    parameter ( p = 2 )
    
    do j = 1, n, 1
      scl = sqrt(sum(x(:,j)**p / (m-1)))
      x(:, j) = x(:,j) / scl
    end do
    
    return
  end subroutine
  
  
  
  ! center and scale simultaneously
  subroutine center_scale(m, n, x)
    implicit none
    ! in/out
    integer :: m, n
    double precision :: x(m, n)
    ! local
    integer :: i, j
    double precision :: mn, var, dt
    
    ! get vectors of column means and sd's
    do j = 1, n, 1
      mn = 0
      var = 0
      do i = 1, m, 1
        dt = x(i,j) - mn
        mn = mn + dt/i
        var = var + dt*(x(i, j) - mn)
      end do
      var = sqrt(var/(m-1))
      
      x(:, j) = (x(:, j) - mn)/var
    end do
    
    return
  end subroutine
  
  
end module
