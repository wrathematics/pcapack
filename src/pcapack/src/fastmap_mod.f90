! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


! purpose
! =======
! this routine computes 
!
!
! notes
! =====
!
! 
!
!
! arguments
! =========
!
!  m, n     (input) integer
!           the number of rows and columns of the input matrix a.
!



! implementation of the fastmap coordinate mapping algorithm:
! "a matrix computation view of fastmap and robustmap dimension reduction algorithms", 
! george ostrouchov, siam j. matrix anal. appl.
!
! x is overwritten with estimates for the first k principal components.
module fastmap_mod
  use :: lapack
  use, intrinsic :: iso_c_binding
  implicit none
  
  
  subroutine cma(n, p, x, k, info)
    !in/out
    integer :: n, p, k, info
    double precision :: x(n, p)
    ! local
    integer :: i, j, l, ncol, allocerr
    double precision, allocatable :: a(:), b(:), y(:), work(:)
    double precision :: tmp, best
    ! parameters
    double precision :: zero, one
    parameter( zero = 0.0d0, one = 1.0d0 )
    
    
    info = 0
    
    ! quick return if possible
    if (k.gt.n) then
      info = -4
      return
    end if
    
    
    allocerr = 0
    
    allocate(a(p), stat=allocerr)
    allocate(b(p), stat=allocerr)
    allocate(y(n), stat=allocerr)
    allocate(work(p), stat=allocerr)
    
    if (allocerr.ne.0) then
      info = 1
      return
    end if
    
    
    ! coordinate mapping algorithm
    do i = 1, k
      ncol = p-i+1
      
      ! select pivot row pair
      call fastmap(n, ncol, i, x(1,i), a, b, work)
      
      ! translate rows to pivot origin
      call dsweep(n, ncol, x(1, i), a, ncol, 2, "-")
      
      ! apply householder reflection on right
      call daxpy(ncol, -1.0d0, a, 1, b, 1)
      
      b(1) = dsgn(b(1)) * dnrm2(ncol, b, 1) + b(1)
      
      tmp = dnrm2(p, b, 1)
      b(1:ncol) = b(1:ncol) / tmp
      
      call dgemv('n', n, ncol, one, x(1,i), n, b, 1, zero, y, 1)
      
      call dger(n, ncol, -2.0d0, y, 1, b, 1, x(1,i), n)
    end do
    
    
    deallocate(a)
    deallocate(b)
    deallocate(y)
    deallocate(work)
    
    
    return
  end
  
  
  
  ! fastmap pivot selection; store the distances
  subroutine fastmap(n, ncol, num, x, a, b, work)
    !in/out
    integer :: n, ncol, num
    double precision :: x(n, ncol), a(ncol), b(ncol), work(ncol)
    ! local
    integer :: i, j, ia, ib
    double precision :: tmp, best
    ! parameters
    double precision :: zero, one, neg1
    parameter( zero = 0.0d0, one = 1.0d0, neg1 = -1.0d0 )
    
    
    ! pick random b in x
    call sample1(1, n, ia)
    
    b = x(ia, :)
    
    
    ! compute all distances in x from b
    best = zero
    
    do i = 1, n
      work = neg1 * x(i, :)
      call daxpy(n, one, b, 1, work, 1)
      tmp = dnrm2(ncol, work, 1)
      
      ! let a in x be the most distant point from b
      if (tmp.gt.best) then
        best = tmp
        ia = i
      end if
    end do
    
    a = x(ia, :)
    
    
    ! compute all distances in x from a
    best = zero
    
    do i = 1, n
      work = neg1 * x(i, :)
      call daxpy(n, one, a, 1, work, 1)
      tmp = dnrm2(ncol, work, 1)
      
      ! let b in x be the most distant point from a
      if (tmp.gt.best) then
        best = tmp
        ia = i
      end if
    end do
    
    b = x(ib, :)
    
    
    return
  end
  
  
  !!!!!$omp parallel if (n > 100)
  !!!!  !$omp do shared(n, m, x) private(j, i) &
  !!!!  !$omp default(none) schedule(static,chunk) &
  !!!!        do i = 1, j-1
  !!!!          do j = 1, n
  !!!!            x(j, i) = x(i, j)
  !!!!          end do
  !!!!        end do 
  !!!!  !$omp end do
  !!!!!$omp end parallel
  
  
  
end module
