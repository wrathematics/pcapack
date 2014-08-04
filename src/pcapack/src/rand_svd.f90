! this source code form is subject to the terms of the mozilla public
! license, v. 2.0. if a copy of the mpl was not distributed with this
! file, you can obtain one at http://mozilla.org/mpl/2.0/.


! copyright 2013, schmidt


! purpose
! =======
! 
! 
!
!
! notes
! =====
!
! the randomized svd method from the paper halko, martinsson, and tropp
! 2011, "finding structure with randomness: probabilistic algorithms 
! for constructing approximate matrix decompositions", siam review 
! 53, 217-288
!
! we implement the method from the paper in 2 ways.  the first is 
! "as is", using a random normal matrix generation for the test
! matrix omega.  additionally, we use the fastmap algorithm to 
! generate the matrix omega.
!
! algorithm is o(m*n*k); see the paper for details.
!
!
! arguments
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


subroutine drandsvd(method, jobu, jobvt, m, n, a, k, q, s, u, vt, info)
  implicit none
  ! in/out
  character*1         method, jobu, jobvt
  integer :: m, n, k, q, info
  double precision :: a(*), s(*), u(*), vt(*)
  ! local
  integer :: i, lwork, rc_min, rc_max
  integer :: allocerr
  double precision :: tmp
  double precision, allocatable, dimension(:) ::  omega, work, y, atmp
  ! parameters
  double precision :: zero, one
  parameter ( zero = 0.0d0, one = 1.0d0 )
  ! subroutines
  intrinsic           min, max
  external            drandsvd_omega_init, drandsvd_stage_a, drandsvd_stage_b
  
  
  info = 0
  allocerr = 0
  
  
  !!! quick return if possible
  if (m < 1 .or. n < 1) then
    return
  end if
  
  if (method /= 'r' .or. method /= 'f') then
    info = -1
    return
  end if
  
  if (jobu /= 'v' .or. jobu /= 'n') then
    info = -2
    return
  end if
  
  if (jobvt /= 'v' .or. jobvt /= 'n') then
    info = -3
    return
  end if
  
  if (k.gt.n .or. k < 1) then
    info = -7
    return
  end if
  
  if (q < 1) then
    info = -9
    return
  end if
  
  rc_min = min(m, n)
  rc_max = max(m, n)
  
  
  !!! allocate and initializeomega=nx2k
  allocate(omega(n*2*k), stat=allocerr)
  call drandsvd_omega_init(method, n, k, omega, info)
  
  
  !!! allocate workspace
!      call dgeqrf(m, n, y, m, tau, tmp, -1, info)
!      lwork = int(tmp)
  ! fixme
  lwork = 10000
  allocate(work(lwork), stat=allocerr)
  
  
  !!! allocate 
  
  allocate(atmp(m*n), stat=allocerr)
  
  ! allocate y=mx2k
  allocate(y(m*2*k), stat=allocerr)
  
  
  !!! stage a
  call drandsvd_stage_a(m, n, k, q, a, y, u, omega, work, lwork, info)
  if (info /= 0) goto 1
  
  
  !!! stage b
  call drandsvd_stage_b(jobu, jobvt, m, n, k, a, y, s, u, vt, work, lwork, info)
  if (info /= 0) goto 1
  
  
  
  !!! free and return
1    continue
  
  write (*,*) info
  
  deallocate(omega)
  deallocate(work)
  deallocate(atmp)
  deallocate(y)
  
  
  return
end




! ------------------------------------------------------
! internal subroutines
! ------------------------------------------------------


! initialize omega, such as via random normal generation.
subroutine drandsvd_omega_init(method, n, k, omega, info)
  implicit none
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
subroutine drandsvd_stage_a(m, n, k, q, a, y, u, omega, work, lwork, info)
  implicit none
  ! in/out
  integer :: m, n, k, q, lwork, info
  double precision :: a(*), y(*), u(*), omega(*), work(*)
  ! local
  integer :: j
  integer :: allocerr
  double precision, allocatable :: tau(:)
  ! parameters
  double precision :: zero, one
  parameter ( zero = 0.0d0, one = 1.0d0 )
  ! subroutines
  external            dgemm, dgesvd, drandsvd_subspaceiter
  
  
  allocate(tau(k), stat=allocerr)
  
  ! y = a * omega
  call dgemm('n', 'n', m, 2*k, n, one, a, m, omega, n, zero, y, m)
  
  ! q = the q part of the qr of y
  call dgeqrf(m, n, y, m, tau, work, lwork, info)
  if (info /= 0) goto 1
  call dorgqr(m, n, k, y, m, tau, work, lwork, info)
  if (info /= 0) goto 1
  
  
  !!! main loop
  
  ! i = 1
  call drandsvd_subspaceiter(m, n, k, q, a, u, y, tau, work, lwork, info)
  if (info /= 0) goto 1
  
  
  ! general
  do j = 2, q
    call drandsvd_subspaceiter(m, n, k, q, a, u, y, tau, work, lwork, info)
    if (info /= 0) goto 1
    
  end do
  
  
  ! free and return
1    continue
  
  deallocate(tau)
  
  return
end



subroutine drandsvd_subspaceiter(m, n, k, q, a, u, y, tau, work, lwork, info)
  implicit none
  ! in/out
  integer :: m, n, k, q, lwork, info
  double precision :: a(*), u(*), y(*), tau(*), work(*)
  ! local
  !
  ! parameters
  double precision :: zero, one
  parameter ( zero = 0.0d0, one = 1.0d0 )
  ! subroutines
  external            dgemm, dgeqrf, dorgqr
  
  
  ! twy_j = a^t * q_j-1
  call dgemm('t', 'n', n, 2*k, n, one, a, m, y, n, zero, y, m)
  
  ! qr of twy_j
  call dgeqrf(m, n, y, m, tau, work, lwork, info)
  if (info /= 0) return
  call dorgqr(m, n, k, y, m, tau, work, lwork, info)
  if (info /= 0) return
  
  ! y = a * q
  call dgemm('n', 'n', m, 2*k, n, one, a, m, y, n, zero, y, m)
  
  ! q = the q part of the qr of y
  call dgeqrf(m, n, y, m, tau, work, lwork, info)
  if (info /= 0) return
  call dorgqr(m, n, k, y, m, tau, work, lwork, info)
  if (info /= 0) return
  
  
  return
end




subroutine drandsvd_stage_b(jobu, jobvt, m, n, k, a, y, s, u, vt, work, lwork, info)
  implicit none
  ! in/out
  character*1         jobu, jobvt
  integer :: m, n, k, info
  double precision :: a(*), y(*), s(*), u(*), vt(*), work(*)
  ! local
  integer :: lwork, qr, qc, br, bc, ur, uc, rc_min
  integer :: allocerr
  double precision, allocatable :: b(:)
  ! parameters
  double precision :: zero, one
  parameter ( zero = 0.0d0, one = 1.0d0 )
  ! subroutines
  external            dgemm, dgesvd
  
  
  rc_min = min(m, n)
  
  br = min(rc_min, 2*k)
  bc = n
  
  qr = m
  qc = bc
  
  ur = m
  uc = min(n, 2*k)
  
  allocate(b(qc * n), stat=allocerr)
  if (allocerr /= 0) goto 1
  
  
  ! b = q^t * a
  call dgemm('t', 'n', qc, n, qr, one, y, qr, a, m, zero, b, qc)
  
  
  ! twu = the u part of the svd of b
  if (jobu == 'v') then
    call dgesvd('o', jobvt, qc, n, b, qc, s, 1.0d0, 1, vt, n, work, lwork, info)
    
    ! u = q * twu
    call dgemm('n', 'n', qr, qc, qc, one, y, qr, b, qc, zero, u, m)
    
    
    ! u = u[, 1:k]
    
    
    if (info /= 0) goto 1
  else
    call dgesvd('n', jobvt, qc, n, b, qc, s, 1.0d0, 1, vt, n,  work, lwork, info)
    
    
    ! vt = vt[1:k, ]
    
    
    if (info /= 0) goto 1
  end if
  
  
  
1    continue
  deallocate(b)
  
  
  return
end




