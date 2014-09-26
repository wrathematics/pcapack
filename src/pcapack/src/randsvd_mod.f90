! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt



! Purpose
! =======
! 
! 
!
!
! Notes
! =====
!
! The randomized svd method from the paper Halko, Martinsson, and Tropp
! 2011, "Finding Structure with Randomness: Probabilistic Algorithms 
! for Constructing Approximate Matrix Decompositions", SIAM Review 
! 53, 217-288.
!
! We implement the method from the paper in 2 ways.  The first is 
! "as is", using a random normal matrix generation for the test
! matrix omega.  Additionally, we use the fastmap algorithm to 
! generate the matrix omega.
!
! Complexity is O(m*n*k); see the paper for details.
!
!
! Arguments
! =========
!
!  m, n     (input) integer
!           the number of rows and columns of the input matrix a.
!
!  a        (input) double precision array
!           input data matrix of size mxn.
!
!  k        (input) integer
!           the desired number of singular vectors.
!
!  u        (output) double precision array
!           mxk
!
!  q        (input) integer
!           
!


! s = min(m, n)
! u = mxk
! vt = 

! omega = nxl, we can take l=2*k for sizes george used
! maybe try using fastmap to create omega




module randsvd_mod
  use :: lapack
  use, intrinsic :: iso_c_binding
  use :: svd_mod
  intrinsic :: min, max
  implicit none
  
  
  integer, public, parameter :: pcapack_randsvd_randinit    = 1
  integer, public, parameter :: pcapack_randsvd_fastmapinit = 2
  
  
  private :: randsvd_omega_init
  public :: randsvd
  
  contains
  
!  subroutine randsvd(nu, nv, m, n, x, s, u, vt, info)
  subroutine randsvd(method, jobu, jobvt, m, n, a, k, q, s, u, vt, info)
    ! in/out
    character*1         method, jobu, jobvt
    integer :: m, n, k, q, info
    double precision :: a(*), s(*), u(*), vt(*)
    ! local
    integer :: i, lwork, rc_min, rc_max
    integer :: allocerr
    double precision :: arrtmp(1), tmp
    double precision, allocatable, dimension(:) ::  omega, work, y, atmp
    
    
    info = 0
    allocerr = 0
    
    
    if (m < 1 .or. n < 1) return
    
    if (method /= 'r' .or. method /= 'f') then
      info = -1
      return
    else if (jobu /= 'v' .or. jobu /= 'n') then
      info = -2
      return
    else if (jobvt /= 'v' .or. jobvt /= 'n') then
      info = -3
      return
    else if (k > n .or. k < 1) then
      info = -7
      return
    else if (q < 1) then
      info = -9
      return
    end if
    
    rc_min = min(m, n)
    rc_max = max(m, n)
    
    
    ! allocate and initializeomega=nx2k
    allocate(omega(n*2*k), stat=allocerr)
    call randsvd_omega_init(method, n, k, omega, info)
    
    
    ! allocate workspace
    call dgeqrf(m, 2*k, y, m, tau, arrtmp, -1, info)
    lwork = int(arrtmp(1))
!    ! fixme
!    lwork = 10000
    allocate(work(lwork), stat=allocerr)
    allocate(atmp(m*n), stat=allocerr)
    allocate(y(m*2*k), stat=allocerr) ! m x (2*k)
    
    
    ! stage a
    call randsvd_stage_a(m, n, k, q, a, y, u, omega, work, lwork, info)
    if (info /= 0) goto 1
    
    
    ! stage b
    call randsvd_stage_b(jobu, jobvt, m, n, k, a, y, s, u, vt, work, lwork, info)
    if (info /= 0) goto 1
    
    
    ! free and return
    1 continue
    write (*,*) info
    
    if (allocated(omega)) deallocate(omega)
    if (allocated(work))  deallocate(work)
    if (allocated(atmp))  deallocate(atmp)
    if (allocated(y))     deallocate(y)
    
    
    return
  end
  
  
  
  
  ! ------------------------------------------------------
  ! internal subroutines
  ! ------------------------------------------------------
  
  
  ! initialize omega, such as via random normal generation.
  subroutine randsvd_omega_init(method, n, k, omega, info)
    ! in/out
    character*1         method
    integer :: n, k, info
    double precision :: omega(*)
    ! subroutines
    external            rnormn
    
    
    ! random (standard) normal initialization
    if (method == 'r') then
      call rnormn(n*2*k, 0.0d0, 1.0d0, omega)
    ! fastmap initialization
    else if (method == 'f') then
      ! fixme
    ! bad argument 'method'
    else
      info = -1
      return
    end if
    
    
    return
  end
  
  
  
  ! stage a from the paper, with the orthonormalization step
  subroutine randsvd_stage_a(m, n, k, q, a, y, u, omega, work, lwork, info)
    ! in/out
    integer :: m, n, k, q, lwork, info
    double precision :: a(*), y(*), u(*), omega(*), work(*)
    ! local
    integer :: j
    integer :: allocerr
    double precision, allocatable :: tau(:)
    
    
    allocate(tau(k), stat=allocerr)
    
    ! y = a * omega
    call dgemm('n', 'n', m, 2*k, n, 1.0d0, a, m, omega, n, 0.0d0, y, m)
    
    ! q = the q part of the qr of y
    call dgeqrf(m, n, y, m, tau, work, lwork, info)
    if (info /= 0) goto 1
    call dorgqr(m, n, k, y, m, tau, work, lwork, info)
    if (info /= 0) goto 1
    
    
    !!! main loop
    
    ! i = 1
    call randsvd_subspaceiter(m, n, k, q, a, u, y, tau, work, lwork, info)
    if (info /= 0) goto 1
    
    
    ! general
    do j = 2, q
      call randsvd_subspaceiter(m, n, k, q, a, u, y, tau, work, lwork, info)
      if (info /= 0) exit
    end do
    
    
    ! free and return
  1    continue
    
    deallocate(tau)
    
    return
  end
  
  
  
  subroutine randsvd_subspaceiter(m, n, k, q, a, u, y, tau, work, lwork, info)
    ! in/out
    integer :: m, n, k, q, lwork, info
    double precision :: a(*), u(*), y(*), tau(*), work(*)
    ! local
    
    
    ! twy_j = a^t * q_j-1
    call dgemm('t', 'n', n, 2*k, n, 1.0d0, a, m, y, n, 0.0d0, y, m)
    
    ! qr of twy_j
    call dgeqrf(m, n, y, m, tau, work, lwork, info)
    if (info /= 0) return
    call dorgqr(m, n, k, y, m, tau, work, lwork, info)
    if (info /= 0) return
    
    ! y = a * q
    call dgemm('n', 'n', m, 2*k, n, 1.0d0, a, m, y, n, 0.0d0, y, m)
    
    ! q = the q part of the qr of y
    call dgeqrf(m, n, y, m, tau, work, lwork, info)
    if (info /= 0) return
    call dorgqr(m, n, k, y, m, tau, work, lwork, info)
    if (info /= 0) return
    
    
    return
  end
  
  
end module

