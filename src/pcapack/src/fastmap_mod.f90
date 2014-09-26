! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


! Purpose
! =======
! cma computes ...
!
!
! Notes
! =====
! Implementation of the fastmap coordinate mapping algorithm:
! "A Matrix Computation View of Fastmap and Robustmap Dimension 
! Reduction Algorithms", George Ostrouchov, SIAM J. Matrix Anal. 
! appl.
! 
!
!
! Arguments
! =========
!
!  m, n     (input) integer
!           the number of rows and columns of the input matrix a.
!
!  x        (in/out) double precision
!           On input, x is the data matrix, and is overwritten with 
!           estimates for the first k principal components.
!
!  k        (input) integer
!           
!
!  info     (output) integer
!           Error numbering.

module fastmap_mod
  use :: lapack
  use :: sgns
  use :: sweeps
  use, intrinsic :: iso_c_binding
  implicit none
  
  
  integer, public, parameter :: pcapack_oom = -2147483647
  
  public :: cma, fastmap
  
  contains
  
  subroutine cma(n, p, x, k, info) &
  bind(C, name='cma_')
    !in/out
    integer, intent(in) :: n, p, k
    integer, intent(out) :: info
    double precision, intent(inout) :: x(n, p)
    ! local
    integer :: i, j, l, ncol, allocerr
    double precision, allocatable :: a(:), b(:), y(:), work(:)
    double precision :: tmp, best
    
    
    info = 0
    
    if (n < 1 .or. p < 1) return
    
    if (k > n) then
      info = -1
      return
    else if (k > p) then
      info = -2
      return
    end if
    
    
    allocerr = 0
    
    allocate(a(p), stat=allocerr)
    allocate(b(p), stat=allocerr)
    allocate(y(n), stat=allocerr)
    allocate(work(p), stat=allocerr)
    
    if (allocerr /= 0) then
      info = pcapack_oom
      goto 1
    end if
    
    
    ! coordinate mapping algorithm
    do i = 1, k
      ncol = p-i+1
      
      ! select pivot row pair
      call fastmap(n, ncol, x(1,i), a, b, work)
      
      ! translate rows to pivot origin
      call sweep(n, ncol, x(1, i), a, ncol, 2, "-")
      
      ! apply householder reflection on right
      call daxpy(ncol, -1.0d0, a, 1, b, 1)
      
      b(1) = sgn(b(1)) * dnrm2(ncol, b, 1) + b(1)
      
      tmp = dnrm2(p, b, 1)
      b(1:ncol) = b(1:ncol) / tmp
      
      call dgemv('n', n, ncol, 1.0d0, x(1,i), n, b, 1, 0.0d0, y, 1)
      
      call dger(n, ncol, -2.0d0, y, 1, b, 1, x(1,i), n)
    end do
    
    
    1 continue
    if (allocated(a)) deallocate(a)
    if (allocated(b)) deallocate(b)
    if (allocated(y)) deallocate(y)
    if (allocated(work)) deallocate(work)
    
    
    return
  end
  
  
  
  ! fastmap pivot selection
    ! Working on the right ncol columns of x
    ! a is the point (in x) furthest from a random vector in x
    ! b is the point (in x) furthest from point a
    ! work is a workspace array
  subroutine fastmap(n, ncol, x, a, b, work) &
  bind(C, name='fastmap_')
    !in/out
    integer, intent(in) :: n, ncol
    double precision, intent(in) :: x(n, ncol)
    double precision, intent(out) :: a(ncol), b(ncol), work(ncol)
    ! local
    integer :: i, j, ia, ib
    double precision :: tmp, best
    
    
    ! pick random b in x
    call sample1(1, n, ia)
    
    b = x(ia, :)
    
    
    ! compute all distances in x from b
    best = 0.0d0
    
    do i = 1, n
      work = -1.0d0 * x(i, :)
      call daxpy(ncol, 1.0d0, b, 1, work, 1)
      tmp = dnrm2(ncol, work, 1)
      
      ! let a in x be the most distant point from b
      if (tmp > best) then
        best = tmp
        ia = i
      end if
    end do
    
    a = x(ia, :)
    
    
    ! compute all distances in x from a
    best = 0.0d0
    
    do i = 1, n
      work = -1.0d0 * x(i, :)
      call daxpy(ncol, 1.0d0, a, 1, work, 1)
      tmp = dnrm2(ncol, work, 1)
      
      ! let b in x be the most distant point from a
      if (tmp > best) then
        best = tmp
        ia = i
      end if
    end do
    
    b = x(ib, :)
    
    
    return
  end
  
end module
